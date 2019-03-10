varying vec4 position;
varying vec3 normal;
varying vec3 light_direction;

void main()
{	
    vec4 ambient = vec4(1, 0, 0, 1);
    vec4 diffuse = vec4(0, 1, 0, 1);
    vec4 specular = vec4(0, 0, 1, 1);
    
    ambient = gl_LightModel.ambient * gl_FrontMaterial.ambient * gl_LightSource[0].ambient;
    float diffuse_clamp = max(0.0, dot(normal, light_direction));
    diffuse = gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse * diffuse_clamp;
    vec3 R = 2.0 * (dot(light_direction, normal) * normal) - light_direction;
    vec4 reflected_direction = vec4(R[0], R[1], R[2], 0);
    float specular_clamp = pow(max(0.0, dot(reflected_direction, -position)), gl_FrontMaterial.shininess);
    specular = gl_LightSource[0].diffuse * gl_FrontMaterial.specular * specular_clamp;
    
    gl_FragColor = ambient + diffuse + specular;
}
