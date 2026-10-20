/**
 * GLSL shaders for point rendering with adaptive sizing
 * (docs/specs/RENDERING.md).
 *
 * Two pairs over one body: the cell cloud carries colour, size and alpha
 * per vertex (every cell differs); a molecule cloud carries position only,
 * with colour, size and alpha as uniforms of its material, so a slider tick
 * changes a uniform and nothing is uploaded.
 *
 * Both implement a smooth-fade distance filter from a world-space anchor
 * (the camera's orbit target / cmd-click recenter point): when
 * uTargetFilterEnabled is on, points fade out over a feather band between
 * uTargetRadius and uTargetRadius + uTargetFeather. Sizes are in CSS
 * pixels: uPixelRatio scales them to the device's.
 */

const vertexBody = `
    uniform float dotSize;
    uniform float uPixelRatio;
    // 1: dotSize is the cells' real size in world units (Visium spots,
    // bins), drawn at that size on screen; uHalfHeight is the canvas's half
    // height in CSS pixels, which turns world units into pixels
    uniform float uWorldSized;
    uniform float uHalfHeight;
    uniform vec3 uTargetCenter;
    uniform float uTargetFilterEnabled;
    varying float vTargetDist;

    void main() {
        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);

        // World-space distance from the orbit target, only when the filter
        // is on: the second matrix multiply is otherwise skipped.
        vTargetDist = 0.0;
        if (uTargetFilterEnabled > 0.5) {
            vec4 worldPos = modelMatrix * vec4(position, 1.0);
            vTargetDist = distance(worldPos.xyz, uTargetCenter);
        }

        // World-space sizing: dotSize is pre-scaled to account for data extent.
        // projectionMatrix[1][1] = 1/tan(fov/2), dividing by -mvPosition.z gives perspective.
        // Points naturally appear bigger when closer and smaller when farther.
        // A real-size cell keeps its size: expression and selection show in
        // its colour and opacity, not by shrinking it
        float px = (uWorldSized > 0.5 ? uHalfHeight : POINT_SCALE) * dotSize
            * projectionMatrix[1][1] / -mvPosition.z;

        // Cull sub-pixel points — without this they jam at the GPU's 1px floor
        // and pile up as overdraw when zoomed far out. Threshold sets the
        // working distance: working distance scales as 1 / threshold, so 0.025
        // gives 40× the camera range of 1.0 before the cloud goes empty. The
        // SM overlay needs the headroom because its molecules sit in [-1,1]
        // while the cell cloud spans [-cs,cs], so the camera distance at the
        // "whole dataset" zoom level is ~cs and SM points hit sub-pixel first.
        if (px < 0.025) {
            gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
            return;
        }
        // Floor at 1 CSS pixel. Many GPUs already round sub-pixel points up
        // to their 1px minimum, but some draw nothing at all, so the default
        // ~0.1px dots left the canvas empty (issue #147). Points below the
        // cull above still drop out; everything else stays visible.
        // A real-size cell may fill much of the view when zoomed in (the GPU
        // still caps a point, typically at 1024 px or more)
        gl_PointSize = clamp(px, 1.0, uWorldSized > 0.5 ? 8192.0 : 200.0) * uPixelRatio;

        gl_Position = projectionMatrix * mvPosition;
    }
`;

const fragmentBody = `
    varying float vTargetDist;
    uniform float uTargetRadius;
    uniform float uTargetFeather;
    uniform float uTargetFilterEnabled;
    uniform float uShape; // 0 = circle (clipped + anti-aliased), 1 = square (full quad)

    void main() {
        // Shape: circles clip outside r=0.5 and smooth the edge; squares
        // render the full point quad.
        float alpha = POINT_ALPHA;
        if (uShape < 0.5) {
            float dist = length(gl_PointCoord - vec2(0.5, 0.5));
            if (dist > 0.5) {
                discard;
            }
            alpha *= smoothstep(0.5, 0.45, dist);
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

        gl_FragColor = vec4(POINT_COLOR, alpha);
    }
`;

/** The cell cloud: colour, size and alpha per vertex. */
export const vertexShader = `
    attribute float size;
    attribute vec3 color;
    attribute float alpha;
    varying vec3 vColor;
    varying float vAlpha;
${vertexBody
  .replace(
    "void main() {",
    "void main() {\n        vColor = color;\n        vAlpha = alpha;",
  )
  .replace("POINT_SCALE", "size")}
`;

export const fragmentShader = `
    varying vec3 vColor;
    varying float vAlpha;
${fragmentBody.replace("POINT_ALPHA", "vAlpha").replace("POINT_COLOR", "vColor")}
`;

/** A molecule cloud: one colour, size and alpha for the whole cloud. */
export const moleculeVertexShader = vertexBody.replace("POINT_SCALE", "1.0");

export const moleculeFragmentShader = `
    uniform vec3 uColor;
    uniform float uAlpha;
${fragmentBody.replace("POINT_ALPHA", "uAlpha").replace("POINT_COLOR", "uColor")}
`;
