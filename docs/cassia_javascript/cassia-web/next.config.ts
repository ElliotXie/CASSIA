import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  experimental: {
    esmExternals: true,
    // Optimize package imports for better tree-shaking
    optimizePackageImports: ['lucide-react', '@radix-ui/react-icons'],
  },
  transpilePackages: [],
  typescript: {
    ignoreBuildErrors: true,
  },
  eslint: {
    ignoreDuringBuilds: true,
  },
  env: {
    NEXT_PUBLIC_SUPABASE_URL: process.env.NEXT_PUBLIC_SUPABASE_URL,
    NEXT_PUBLIC_SUPABASE_ANON_KEY: process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY,
  },
  // Tree-shake icon libraries
  modularizeImports: {
    'lucide-react': {
      transform: 'lucide-react/dist/esm/icons/{{member}}',
    },
  },
  // Production optimizations
  compiler: {
    // Remove console logs in production (keep error and warn)
    removeConsole: process.env.NODE_ENV === 'production' ? {
      exclude: ['error', 'warn'],
    } : false,
  },
  // Enable compression
  compress: true,
};

// Bundle analyzer - run with: ANALYZE=true npm run build
const withBundleAnalyzer = require('@next/bundle-analyzer')({
  enabled: process.env.ANALYZE === 'true',
});

export default withBundleAnalyzer(nextConfig);
