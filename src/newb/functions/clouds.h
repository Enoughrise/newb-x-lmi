#ifndef CLOUDS_H
#define CLOUDS_H

#include "noise.h"
#include "sky.h"

// simple clouds 2D noise
float cloudNoise2D(vec2 p, highp float t, float rain) {
  t *= NL_CLOUD1_SPEED;
  p += t;
  p.y += 3.0*sin(0.3*p.x + 0.1*t);

  vec2 p0 = floor(p);
  vec2 u = p-p0;
  u *= u*(3.0-2.0*u);
  vec2 v = 1.0-u;

  // rain transition
  float n = mix(
    mix(rand(p0),rand(p0+vec2(1.0,0.0)), u.x),
    mix(rand(p0+vec2(0.0,1.0)),rand(p0+vec2(1.0,1.0)), u.x),
    u.y
  );
  n *= 0.5 + 0.5*sin(p.x*0.6 - 0.5*t)*sin(p.y*0.6 + 0.8*t);
  n = min(n*(1.0+rain), 1.0);
  return n*n;
}

// simple clouds
vec4 renderCloudsSimple(nl_skycolor skycol, vec3 pos, highp float t, float rain) {
  pos.xz *= NL_CLOUD1_SCALE;

  float d = cloudNoise2D(pos.xz, t, rain);
  vec4 col = vec4(skycol.horizonEdge + skycol.zenith, smoothstep(0.1,0.6,d));
  col.rgb += 1.5*dot(col.rgb, vec3(0.3,0.4,0.3))*smoothstep(0.6,0.2,d)*col.a;
  col.rgb *= 1.0 - 0.8*rain;
  return col;
}

// rounded clouds

// rounded clouds 3D density map
#if NL_CLOUD2_TYPE == 1
float cloudDf(vec3 pos, float rain, float time, vec2 boxiness) {
  boxiness *= 0.999;
  vec2 p0 = floor(pos.xz);
  vec2 u = max((pos.xz-p0-boxiness)/(1.0-boxiness), 0.0);
  u *= u*(3.0 - 2.0*u);
  //vec2 v = 1.0 - u;
#elif NL_CLOUD2_TYPE == 2
// 2d noise
float noise(vec2 p){
  vec2 p0 = floor(p);
	vec2 u = p-p0;

	u *= u*(3.0-2.0*u);
	vec2 v = 1.0 - u;

	float c1 = rand(p0);
	float c2 = rand(p0+vec2(1.0,0.0));
	float c3 = rand(p0+vec2(0.0,1.0));
	float c4 = rand(p0+vec2(1.0,1.0));

	float n = v.y*mix(c1,c2,u.x) + u.y*(c3*v.x+c4*u.x);
	return n;
}
float cloudDf(vec3 pos, float rain, float time, vec2 boxiness) {
  pos.x += 0.35*noise(6.35*pos.xz) + 0.05*time;
  pos.y += 0.0*noise(0.0*pos.xz);
  pos.z += 0.35*noise(6.35*pos.xz) + 0.05*time;
  boxiness *= 0.999;
  vec2 p0 = floor(pos.xz);
  vec2 u = max((pos.xz-p0-boxiness)/(1.0-boxiness), 0.0);
  //u *= u*(3.0 - 2.0*u);
  vec2 v = 1.0 - u;
#endif

  vec4 r = vec4(rand(p0), rand(p0+vec2(1.0,0.0)), rand(p0+vec2(1.0,1.0)), rand(p0+vec2(0.0,1.0)));
  r = smoothstep(0.1001+0.2*rain, 0.1+0.2*rain*rain, r); // rain transition

  float n = mix(mix(r.x,r.y,u.x), mix(r.w,r.z,u.x), u.y);

  // round y
  n *= 1.0 - 1.5*smoothstep(boxiness.y, 2.0 - boxiness.y, 2.0*abs(pos.y-0.55));
  n = max(1.0*(n-0.0), 0.0); // smoothstep(0.2, 1.0, n)
  n *= n*(3.0 - 2.0*n);
  return n;
}

vec4 renderClouds(
    vec3 vDir, vec3 vPos, float rain, float time, vec3 horizonCol, vec3 zenithCol,
    const int steps, const float thickness, const float thickness_rain, const float speed,
    const vec2 scale, const float density, const vec2 boxiness
) {
  float height = 7.0*mix(thickness, thickness_rain, rain);
  float stepsf = float(steps);

  // scaled ray offset
  vec3 deltaP;
  deltaP.y = 1.0;
  deltaP.xz = height*scale*vDir.xz/(0.02+0.98*abs(vDir.y));

  // local cloud pos
  vec3 pos;
  pos.y = 0.0;
  pos.xz = scale*(vPos.xz + vec2(1.0,0.5)*(time*speed));
  pos += deltaP;

  deltaP /= -stepsf;

  // alpha, gradient
  vec2 d = vec2(0.0,1.0);
  for (int i=1; i<=steps; i++) {
    float m = cloudDf(pos, rain, time, boxiness);
    d.x += m;
    d.y = mix(d.y, pos.y, m);
    pos += deltaP;
  }
  
  d.x *= smoothstep(0.7, 1.4, d.x);
  d.x /= (stepsf/density) + d.x;

  if (vPos.y < 0.0) { // view from top
    d.y = 1.0 - d.y;
  }

  vec4 col = vec4(horizonCol + zenithCol, d.x);
  col.rgb *= mix(0.28,1.01,d.y);
  col.rgb += dot(col.rgb, vec3(0.3,0.4,0.3))*d.y*d.y;
  col.rgb *= 1.0 - 0.8*rain;

  return col;
}

// aurora is rendered on clouds layer
#ifdef NL_AURORA
vec4 renderAurora(vec3 p, float t, float rain, vec3 FOG_COLOR) {
  t *= NL_AURORA_VELOCITY;
  p.xz *= NL_AURORA_SCALE;
  p.xz += 0.05*sin(p.x*4.0 + 20.0*t);

  float d0 = sin(p.x*0.25 + t + sin(p.z*0.9));
  float d1 = sin(p.z*0.1 - t + sin(p.x*0.2));
  float d2 = sin(p.z*0.25 + 1.0*sin(d0 + d1*2.0) + d1*2.0 + d0*1.0);
  d0 *= d0; d1 *= d1; d2 *= d2;
  d2 = d0/(1.0 + d2/NL_AURORA_WIDTH);

  float mask = (1.0-0.8*rain)*max(1.0 - 4.0*max(FOG_COLOR.b, FOG_COLOR.g), 0.0);
  return vec4(NL_AURORA*mix(NL_AURORA_COL1,NL_AURORA_COL2,d1),1.0)*d2*mask;
}
#endif

#endif
