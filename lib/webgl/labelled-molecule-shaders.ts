/**
 * Shaders for the labelled single-molecule viewer.
 *
 * Every molecule carries one category index per column (gene / domain / cell)
 * as a static vertex attribute. Selection state lives in small lookup textures
 * — one texel per category — so toggling a checkbox re-uploads well under a
 * kilobyte instead of rebuilding colour/size/alpha arrays for every point.
 *
 * A molecule is drawn "selected" when it passes ALL three menus. Each LUT holds
 * 1.0 for a checked value, and is filled entirely with 1.0 when its menu is
 * empty — so an empty menu contributes no constraint and the product still
 * reads 1.0.
 *
 * The palette LUT packs colour and size together: RGB is the value's colour and
 * A is its size multiplier, encoded over [0, SIZE_ENCODE_RANGE]. Packing them
 * into one RGBA8 texture avoids depending on float-texture support.
 */

/** Per-value size multipliers are encoded into an 8-bit alpha over this range. */
export const SIZE_ENCODE_RANGE = 4.0;

export const labelledMoleculeVertexShader = `
    attribute float aGene;
    attribute float aDomain;
    attribute float aCell;
    attribute float aSize;
    attribute float aAlpha;

    uniform sampler2D uSelGene;
    uniform sampler2D uSelDomain;
    uniform sampler2D uSelCell;
    uniform sampler2D uPalette;

    uniform float uNGene;
    uniform float uNDomain;
    uniform float uNCell;
    uniform float uNColorBy;
    uniform float uColorBy;          // 0 = gene, 1 = domain, 2 = cell

    uniform float dotSize;
    uniform float uGlobalSize;
    uniform float uGlobalAlpha;
    uniform float uUnselectedAlpha;
    uniform float uUnselectedHidden; // 1.0 = cull unselected entirely
    uniform float uUnselectedSize;   // size multiplier for greyed points
    uniform vec3 uUnselectedColor;

    varying vec3 vColor;
    varying float vAlpha;

    // One texel of a width-n 1D lookup texture, sampled at its centre.
    float lut1(sampler2D tex, float index, float n) {
        return texture2D(tex, vec2((index + 0.5) / n, 0.5)).r;
    }

    void main() {
        float selected = lut1(uSelGene, aGene, uNGene)
                       * lut1(uSelDomain, aDomain, uNDomain)
                       * lut1(uSelCell, aCell, uNCell);

        if (selected < 0.5 && uUnselectedHidden > 0.5) {
            gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
            return;
        }

        // Index the palette by whichever column is currently colouring.
        float colorIndex = uColorBy < 0.5 ? aGene
                         : (uColorBy < 1.5 ? aDomain : aCell);
        vec4 entry = texture2D(uPalette,
                     vec2((colorIndex + 0.5) / uNColorBy, 0.5));

        vColor = selected > 0.5 ? entry.rgb : uUnselectedColor;
        vAlpha = aAlpha * uGlobalAlpha
               * (selected > 0.5 ? 1.0 : uUnselectedAlpha);

        // Per-value size only applies to selected points; greyed ones fall back
        // to a uniform multiplier so a big category can't dominate the backdrop.
        float valueSize = entry.a * ${SIZE_ENCODE_RANGE.toFixed(1)};
        float sizeMul = selected > 0.5 ? valueSize : uUnselectedSize;

        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);

        gl_PointSize = aSize * sizeMul * uGlobalSize * dotSize
                     * projectionMatrix[1][1] / -mvPosition.z;

        // Cull sub-pixel points; without this they jam at the GPU's 1px floor
        // and pile up as overdraw when zoomed out (same threshold as the
        // single-cell shader).
        if (gl_PointSize < 0.025) {
            gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
            return;
        }
        gl_PointSize = min(gl_PointSize, 200.0);

        gl_Position = projectionMatrix * mvPosition;
    }
`;

export const labelledMoleculeFragmentShader = `
    uniform float uShape; // 0 = circle, 1 = square

    varying vec3 vColor;
    varying float vAlpha;

    void main() {
        float alpha = vAlpha;

        if (uShape < 0.5) {
            float dist = length(gl_PointCoord - vec2(0.5));

            if (dist > 0.5) discard;
            alpha *= smoothstep(0.5, 0.45, dist);
        }

        if (alpha < 0.004) discard;

        gl_FragColor = vec4(vColor, alpha);
    }
`;
