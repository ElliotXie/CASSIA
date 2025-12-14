import createNextIntlPlugin from 'next-intl/plugin';

const withNextIntl = createNextIntlPlugin('./i18n/request.ts');

/** @type {import('next').NextConfig} */
const nextConfig = {
  eslint: {
    ignoreDuringBuilds: true,
  },
  typescript: {
    ignoreBuildErrors: true,
  },
  images: {
    formats: ['image/avif', 'image/webp'],
  },
  async headers() {
    return [
      {
        source: '/images/:path*',
        headers: [
          {
            key: 'Cache-Control',
            value: 'public, max-age=31536000, immutable',
          },
        ],
      },
      {
        source: '/_next/static/:path*',
        headers: [
          {
            key: 'Cache-Control',
            value: 'public, max-age=31536000, immutable',
          },
        ],
      },
    ]
  },
  async redirects() {
    return [
      // Redirect old docs URLs to new R version
      {
        source: '/:locale/docs/:slug',
        destination: '/:locale/docs/r/:slug',
        permanent: true,
      },
      // Redirect old vignette URLs to new R version
      {
        source: '/:locale/vignette/:slug',
        destination: '/:locale/vignette/r/:slug',
        permanent: true,
      },
    ]
  },
}

export default withNextIntl(nextConfig)
