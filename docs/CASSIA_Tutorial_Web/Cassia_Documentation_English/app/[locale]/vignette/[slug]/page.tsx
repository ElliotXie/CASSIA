import { getVignetteBySlug, getAllVignettes } from "@/lib/vignette"
import DocPageClient from "../../docs/[slug]/DocPageClient"
import { notFound } from "next/navigation"
import { routing } from "@/i18n/routing"

interface VignettePageProps {
  params: Promise<{
    slug: string
    locale: string
  }>
}

// This is the Server Component for Vignette pages
export default async function VignettePage({ params }: VignettePageProps) {
  const { slug, locale } = await params

  // Fetch data on the server with locale
  const vignette = getVignetteBySlug(slug, locale)

  // Handle case where vignette is not found
  if (!vignette) {
    notFound()
  }

  // Ensure the vignette data structure matches what DocPageClient expects
  const processedVignette = {
    slug: vignette.slug,
    frontmatter: {
      ...vignette.frontmatter,
      // Ensure title property exists (required by DocType)
      title: vignette.frontmatter.title || "Untitled Vignette"
    },
    content: vignette.content,
    isFallback: vignette.isFallback
  }

  // Pass the fetched data and params to the DocPageClient
  return <DocPageClient doc={processedVignette} params={{ slug, locale }} />
}

// Add this function to generate static paths at build time
export async function generateStaticParams() {
  const params: { locale: string; slug: string }[] = []

  for (const locale of routing.locales) {
    const vignettes = getAllVignettes(locale)
    for (const vignette of vignettes) {
      params.push({
        locale,
        slug: vignette.slug,
      })
    }
  }

  return params
}
