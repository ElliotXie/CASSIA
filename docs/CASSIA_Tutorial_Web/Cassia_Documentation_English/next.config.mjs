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
    unoptimized: true,
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
