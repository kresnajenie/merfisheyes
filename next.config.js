/** @type {import('next').NextConfig} */
const nextConfig = {
  eslint: {
    ignoreDuringBuilds: true,
  },
  experiments: {
    asyncWebAssembly: true,
    // other experiments if needed
  },
};

module.exports = nextConfig;
