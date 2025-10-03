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

        // Calculate distance from camera
        float distance = -mvPosition.z;
        vDistance = distance;

        // Get base size from the size attribute, scaled by the dotSize uniform
        float baseSize = size * dotSize * 0.4; // Scale factor to make it reasonable

        // Dynamic sizing based on distance with smoother transitions
        float minSize = max(0.5, dotSize * 0.2); // Minimum size scales with dotSize
        float maxSize = min(50.0, dotSize * 6.0); // Maximum size scales with dotSize
        float zoomFactor = 150.0; // LOWER value makes points shrink faster when zooming in

        // Use a smooth curve for size transition based on distance
        // This creates a more natural zoom feeling
        float distanceRatio = zoomFactor / distance;

        // Smooth adaptive sizing with cubic easing
        float t = clamp((distance - 100.0) / 200.0, 0.0, 1.0); // Shorter distance range for faster transition
        float easedT = 1.0 - (1.0 - t) * (1.0 - t) * (1.0 - t); // Cubic ease-out

        // Blend between close-up and far-away behaviors
        float closeUpFactor = 1.0;  // Size multiplier when close to camera
        float farAwayFactor = 2.0;   // Size multiplier when far from camera
        float scaleFactor = mix(closeUpFactor, farAwayFactor, easedT);

        // Calculate final adaptive size
        float adaptiveSize = baseSize * distanceRatio * scaleFactor;

        // Clamp size between min and max
        gl_PointSize = clamp(adaptiveSize, minSize, maxSize);
        gl_Position = projectionMatrix * mvPosition;
    }
`;

export const fragmentShader = `
    varying vec3 vColor;
    varying float vAlpha;
    varying float vDistance;

    void main() {
        // Create circular points instead of squares
        float dist = length(gl_PointCoord - vec2(0.5, 0.5));
        if (dist > 0.5) {
            discard;
        }

        // Enhanced edge effect for all points
        float edgeWidth = 0.15;  // Wider edge
        float distFromCenter = dist;

        // Smooth edge effect that transitions based on distance
        float edgeEffect = 1.0;
        float edgeFactor = smoothstep(0.5 - edgeWidth, 0.5, distFromCenter);

        // Transition edge effect based on distance
        float distanceFactor = smoothstep(150.0, 50.0, vDistance);
        edgeEffect = mix(1.0, 0.7, edgeFactor * distanceFactor);

        // Add subtle anti-aliasing at the edge
        float alpha = vAlpha;
        if (dist > 0.48) {
            alpha *= smoothstep(0.5, 0.48, dist);
        }

        // Apply edge effect to color
        vec3 finalColor = mix(vec3(1.0, 1.0, 1.0), vColor, edgeEffect);

        gl_FragColor = vec4(finalColor, alpha);
    }
`;
