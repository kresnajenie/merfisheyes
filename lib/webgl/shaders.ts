/**
 * Custom GLSL shaders for point cloud rendering with adaptive sizing
 */

export const vertexShader = `
    attribute float size;
    attribute vec3 color;
    attribute float alpha;
    uniform float dotSize;
    varying vec3 vColor;
    varying float vAlpha;
    varying float vDistance;

    void main() {
        vColor = color;
        vAlpha = alpha;
        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);

        // Distance from camera (for fragment shader edge effects)
        vDistance = -mvPosition.z;

        // World-space sizing: dotSize is pre-scaled to account for data extent.
        // projectionMatrix[1][1] = 1/tan(fov/2), dividing by -mvPosition.z gives perspective.
        // Points naturally appear bigger when closer and smaller when farther.
        // size attribute = per-point multiplier (expression/slider).
        gl_PointSize = size * dotSize * projectionMatrix[1][1] / -mvPosition.z;
        gl_PointSize = clamp(gl_PointSize, 0.5, 200.0);

        gl_Position = projectionMatrix * mvPosition;
    }
`;

export const fragmentShader = `
    varying vec3 vColor;
    varying float vAlpha;
    varying float vDistance;

    void main() {
        // Circular points with smooth anti-aliased edges
        float dist = length(gl_PointCoord - vec2(0.5, 0.5));
        if (dist > 0.5) {
            discard;
        }

        // Anti-alias the edge
        float alpha = vAlpha * smoothstep(0.5, 0.45, dist);

        gl_FragColor = vec4(vColor, alpha);
    }
`;
