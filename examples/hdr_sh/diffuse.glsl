// GLSL 330-compatible helper. Upload coefficient-major canonical RGB rows.
// Coefficients already include cosine convolution exactly once.
// n must be a unit normal in the same world frame as the rotated environment.
vec3 evaluateSHIrradiance(vec3 n, vec3 irradianceSH[9]) {
    float x = n.x, y = n.y, z = n.z;
    return irradianceSH[0] * 0.28209479177387814
         - irradianceSH[1] * (0.4886025119029199 * y)
         + irradianceSH[2] * (0.4886025119029199 * z)
         - irradianceSH[3] * (0.4886025119029199 * x)
         + irradianceSH[4] * (1.0925484305920792 * x * y)
         - irradianceSH[5] * (1.0925484305920792 * y * z)
         + irradianceSH[6] * (0.31539156525252005 * (3.0 * z * z - 1.0))
         - irradianceSH[7] * (1.0925484305920792 * x * z)
         + irradianceSH[8] * (0.5462742152960396 * (x * x - y * y));
}

vec3 evaluateSHLambertian(vec3 n, vec3 albedoLinear, vec3 irradianceSH[9]) {
    return albedoLinear * evaluateSHIrradiance(n, irradianceSH)
         * 0.3183098861837907; // 1/pi; no additional convolution
}

// Use once at final output only, and only without framebuffer sRGB conversion.
vec3 linearToSRGB(vec3 linearRGB) {
    vec3 low = 12.92 * linearRGB;
    vec3 high = 1.055 * pow(max(linearRGB, vec3(0.0)), vec3(1.0 / 2.4)) - 0.055;
    return mix(high, low, lessThanEqual(linearRGB, vec3(0.0031308)));
}

// Match the headless sphere PNGs: linear exposure, Reinhard, then sRGB.
vec3 referenceDisplayRGB(vec3 lambertianLinear, float exposure) {
    vec3 exposed = max(lambertianLinear, vec3(0.0)) * exposure;
    return linearToSRGB(exposed / (vec3(1.0) + exposed));
}
