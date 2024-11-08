$input v_texcoord0

#include <bgfx_shader.sh>
#include <newb/functions/tonemap.h>

#ifndef INSTANCING
  #include <newb/config.h>

  uniform vec4 SunMoonColor;

  SAMPLER2D_AUTOREG(s_SunMoonTexture);
#endif

void main() {
  #ifndef INSTANCING
    vec4 color = texture2D(s_SunMoonTexture, v_texcoord0);
    color.rgb *= SunMoonColor.rgb;
    color.rgb *= 4.4*color.rgb;
    float tr = 1.0 - SunMoonColor.a;
    color.a *= 1.0 - tr*tr*tr*tr*tr;
    color.rgb = colorCorrection(color.rgb);
    gl_FragColor = color;
  #else
    gl_FragColor = vec4(0.0, 0.0, 0.0, 0.0);
  #endif
}
