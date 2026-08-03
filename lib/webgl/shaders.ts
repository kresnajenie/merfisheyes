/**
 * Custom GLSL shaders for point cloud rendering with adaptive sizing.
 *
 * Also implements a smooth-fade distance filter from a world-space anchor
 * (the camera's orbit target / cmd-click recenter point). When
 * uTargetFilterEnabled is on, points fade out over a feather band between
 * uTargetRadius and uTargetRadius + uTargetFeather.
 */

export const vertexShader = `
    attribute float size;
    attribute vec3 color;
    attribute float alpha;
    uniform float dotSize;
    uniform vec3 uTargetCenter;
    varying vec3 vColor;
    varying float vAlpha;
    varying float vDistance;
    varying float vTargetDist;

    void main() {
        vColor = color;
        vAlpha = alpha;
        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);

        // Distance from camera (for fragment shader edge effects)
        vDistance = -mvPosition.z;

        // World-space distance from the orbit target (for the distance filter).
        vec4 worldPos = modelMatrix * vec4(position, 1.0);
        vTargetDist = distance(worldPos.xyz, uTargetCenter);

        // World-space sizing: dotSize is pre-scaled to account for data extent.
        // projectionMatrix[1][1] = 1/tan(fov/2), dividing by -mvPosition.z gives perspective.
        // Points naturally appear bigger when closer and smaller when farther.
        // size attribute = per-point multiplier (expression/slider).
        gl_PointSize = size * dotSize * projectionMatrix[1][1] / -mvPosition.z;

        // Cull sub-pixel points — without this they jam at the GPU's 1px floor
        // and pile up as overdraw when zoomed far out. Threshold sets the
        // working distance: working distance scales as 1 / threshold, so 0.025
        // gives 40× the camera range of 1.0 before the cloud goes empty. The
        // SM overlay needs the headroom because its molecules sit in [-1,1]
        // while the cell cloud spans [-cs,cs], so the camera distance at the
        // "whole dataset" zoom level is ~cs and SM points hit sub-pixel first.
        if (gl_PointSize < 0.025) {
            gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
            return;
        }
        gl_PointSize = min(gl_PointSize, 200.0);

        gl_Position = projectionMatrix * mvPosition;
    }
`;

export const fragmentShader = `
    varying vec3 vColor;
    varying float vAlpha;
    varying float vDistance;
    varying float vTargetDist;
    uniform float uTargetRadius;
    uniform float uTargetFeather;
    uniform float uTargetFilterEnabled;
    uniform float uShape; // 0 = circle (clipped + anti-aliased), 1 = square (full quad)

    void main() {
        // Shape: circles clip outside r=0.5 and smooth the edge; squares
        // render the full point quad (alphaTest handles edges in opaque mode).
        float dist = length(gl_PointCoord - vec2(0.5, 0.5));
        float alpha;
        if (uShape < 0.5) {
            if (dist > 0.5) {
                discard;
            }
            alpha = vAlpha * smoothstep(0.5, 0.45, dist);
        } else {
            alpha = vAlpha;
        }

        // Distance-from-target fade: full alpha inside uTargetRadius,
        // smoothly to 0 across the feather band.
        if (uTargetFilterEnabled > 0.5) {
            float fade = 1.0 - smoothstep(
                uTargetRadius,
                uTargetRadius + uTargetFeather,
                vTargetDist
            );
            alpha *= fade;
            if (alpha < 0.005) discard;
        }

        gl_FragColor = vec4(vColor, alpha);
    }
`;
