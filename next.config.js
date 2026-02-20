/** @type {import('next').NextConfig} */
const nextConfig = {
  eslint: {
    ignoreDuringBuilds: true,
  },
  experiments: {
    asyncWebAssembly: true,
    // other experiments if needed
  },
  images: {
    remotePatterns: [
      { hostname: "lh3.googleusercontent.com" },
    ],
  },
};

module.exports = nextConfig;
