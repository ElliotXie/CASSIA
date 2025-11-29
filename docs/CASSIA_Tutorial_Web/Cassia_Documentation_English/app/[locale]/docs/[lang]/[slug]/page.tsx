import { getDocBySlug, getAllDocs, type ProgrammingLang } from "@/lib/docs"
import DocPageClient from "./DocPageClient"
import { notFound } from "next/navigation"
import { routing } from "@/i18n/routing"

interface DocPageProps {
  params: Promise<{
    slug: string
    locale: string
    lang: ProgrammingLang
  }>
}

// This is the Server Component for Doc pages
export default async function DocPage({ params }: DocPageProps) {
  const { slug, locale, lang } = await params

  // Validate programming language
  if (lang !== 'r' && lang !== 'python') {
    notFound()
  }

  // Fetch data on the server with locale and programming language
  const doc = getDocBySlug(slug, locale, lang)

  // Handle case where doc is not found
  if (!doc) {
    notFound()
  }

  // Ensure the doc data structure has the required title
  const processedDoc = {
    slug: doc.slug,
    frontmatter: {
      ...doc.frontmatter,
      // Ensure title property exists (required by DocType)
      title: doc.frontmatter.title || "Untitled Document"
    },
    content: doc.content,
    isFallback: doc.isFallback,
    isPythonFallback: doc.isPythonFallback
  }

  // Pass the fetched data and params to the DocPageClient
  return <DocPageClient doc={processedDoc} params={{ slug, locale, lang }} />
}

// Add this function to generate static paths at build time
export async function generateStaticParams() {
  const params: { locale: string; lang: ProgrammingLang; slug: string }[] = []
  const programmingLanguages: ProgrammingLang[] = ['r', 'python']

  for (const locale of routing.locales) {
    for (const lang of programmingLanguages) {
      const docs = getAllDocs(locale, lang)
      for (const doc of docs) {
        params.push({
          locale,
          lang,
          slug: doc.slug,
        })
      }
    }
  }

  return params
}
