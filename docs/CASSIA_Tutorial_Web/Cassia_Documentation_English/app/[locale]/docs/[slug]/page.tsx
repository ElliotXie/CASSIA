import { getDocBySlug, getAllDocs } from "@/lib/docs"
import DocPageClient from "./DocPageClient"
import { notFound } from "next/navigation"
import { routing } from "@/i18n/routing"

interface DocPageProps {
  params: Promise<{
    slug: string
    locale: string
  }>
}

// This is the Server Component for Doc pages
export default async function DocPage({ params }: DocPageProps) {
  const { slug, locale } = await params

  // Fetch data on the server with locale
  const doc = getDocBySlug(slug, locale)

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
    isFallback: doc.isFallback
  }

  // Pass the fetched data and params to the DocPageClient
  return <DocPageClient doc={processedDoc} params={{ slug, locale }} />
}

// Add this function to generate static paths at build time
export async function generateStaticParams() {
  const params: { locale: string; slug: string }[] = []

  for (const locale of routing.locales) {
    const docs = getAllDocs(locale)
    for (const doc of docs) {
      params.push({
        locale,
        slug: doc.slug,
      })
    }
  }

  return params
}
