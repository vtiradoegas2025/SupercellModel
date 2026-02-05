// Common lighting functions for shaders

// Simple Phong lighting
vec3 phongLighting(vec3 normal, vec3 viewDir, vec3 lightDir, vec3 color) 
{
    vec3 ambient = color * 0.3;
    
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = color * diff * 0.7;
    
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);
    vec3 specular = vec3(1.0) * spec * 0.2;
    
    return ambient + diffuse + specular;
}

// Calculate normal from gradient (for volume rendering)
vec3 calculateNormal(sampler3D volume, vec3 pos, float stepSize) 
{
    vec3 eps = vec3(stepSize, 0.0, 0.0);
    float dx = texture(volume, pos + eps.xyz).a - texture(volume, pos - eps.xyz).a;
    float dy = texture(volume, pos + eps.yxz).a - texture(volume, pos - eps.yxz).a;
    float dz = texture(volume, pos + eps.yzx).a - texture(volume, pos - eps.yzx).a;
    
    vec3 normal = normalize(vec3(dx, dy, dz));
    return normal;
}
