import { getVignetteBySlug, getAllVignettes, type ProgrammingLang } from "@/lib/vignette"
import DocPageClient from "../../../docs/[lang]/[slug]/DocPageClient"
import { notFound } from "next/navigation"
import { routing } from "@/i18n/routing"

interface VignettePageProps {
  params: Promise<{
    slug: string
    locale: string
    lang: ProgrammingLang
  }>
}

// This is the Server Component for Vignette pages
export default async function VignettePage({ params }: VignettePageProps) {
  const { slug, locale, lang } = await params

  // Validate programming language
  if (lang !== 'r' && lang !== 'python') {
    notFound()
  }

  // Fetch data on the server with locale and programming language
  const vignette = getVignetteBySlug(slug, locale, lang)

  // Handle case where vignette is not found
  if (!vignette) {
    notFound()
  }

  // Ensure the vignette data structure matches what DocPageClient expects
  const processedVignette = {
    slug: vignette.slug,
    frontmatter: {
      ...vignette.frontmatter,
      title: vignette.frontmatter.title || "Untitled Vignette"
    },
    content: vignette.content,
    isFallback: vignette.isFallback,
    isPythonFallback: vignette.isPythonFallback
  }

  // Pass the fetched data and params to the DocPageClient
  return <DocPageClient doc={processedVignette} params={{ slug, locale, lang }} />
}

// Add this function to generate static paths at build time
export async function generateStaticParams() {
  const params: { locale: string; lang: ProgrammingLang; slug: string }[] = []
  const programmingLanguages: ProgrammingLang[] = ['r', 'python']

  for (const locale of routing.locales) {
    for (const lang of programmingLanguages) {
      const vignettes = getAllVignettes(locale, lang)
      for (const vignette of vignettes) {
        params.push({
          locale,
          lang,
          slug: vignette.slug,
        })
      }
    }
  }

  return params
}
