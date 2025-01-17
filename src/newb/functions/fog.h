#ifndef FOG_H
#define FOG_H

#include "detection.h"
#include "sky.h"

float nlRenderFogFade(nl_skycolor skycol, float relativeDist, nl_environment env, vec3 FOG_COLOR, vec2 FOG_CONTROL, vec3 wPos, vec3 tsp, float time) {
  #ifdef NL_FOG
    float fade = smoothstep(FOG_CONTROL.x, FOG_CONTROL.y, relativeDist);

    // misty effect
    float density = NL_MIST_DENSITY*(19.0 - 18.0*FOG_COLOR.g);
    fade += (1.0-fade)*(0.3-0.3*exp(-relativeDist*relativeDist*density));

    return NL_FOG * fade;
  #else
    return 0.0;
  #endif
}

float nlRenderGodRay(vec3 cPos, vec3 worldPos, float t, vec2 uv1, float relativeDist, vec3 FOG_COLOR, float fogColor) {
  vec3 offset = cPos - 16.0*fract(worldPos*0.0625);
  offset = abs(2.0*fract(offset*0.0625)-1.0);
  offset = offset*offset*(3.0-2.0*offset);

  vec3 nrmof = normalize(worldPos);

  float u = nrmof.z/length(nrmof.zy);
  float diff = dot(offset,vec3(0.1,0.2,1.0)) + 0.07*t;
  float mask = nrmof.x*nrmof.x;

  float vol = sin(7.0*u + 1.5*diff)*sin(3.0*u + diff);
  vol += sin(5.0*u + 0.4*diff)*sin(4.0*u + 0.7*diff);
  vol *= vol*mask*uv1.y;
  vol *= min(7.0*relativeDist*(1.0-mask),1.0);
  vol *= clamp(3.0*(FOG_COLOR.r-FOG_COLOR.b),0.0,1.0);
  vol = clamp(vol,0.0,1.0);
  
  return fogColor + 0.3*vol*(NL_GODRAY-fogColor);

  vol *= clamp(3.0*(FOG_COLOR.r-FOG_COLOR.b), 0.0, 1.0);

  vol = smoothstep(0.0, 0.1, vol);
  return vol;
}

#endif
