import { getDocBySlug, getAllDocs } from "@/lib/docs"
import DocPageClient from "./DocPageClient"
import { notFound } from "next/navigation"

interface DocPageProps {
  params: {
    slug: string
  }
}

// This is the Server Component for Doc pages
export default async function DocPage({ params }: DocPageProps) {
  // Fetch data on the server
  const doc = getDocBySlug(params.slug)

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
    content: doc.content
  }

  // Pass the fetched data and params to the DocPageClient
  return <DocPageClient doc={processedDoc} params={params} />
}

// Add this function to generate static paths at build time
export async function generateStaticParams() {
  const docs = getAllDocs()
  return docs.map((doc) => ({
    slug: doc.slug,
  }))
}

// Optional: Add generateMetadata to set page title, etc.
// import { Metadata } from 'next'
// export async function generateMetadata({ params }: DocPageProps): Promise<Metadata> {
//   const doc = getDocBySlug(params.slug)
//   if (!doc) {
//     return { title: "Not Found" }
//   }
//   return {
//     title: doc.frontmatter.title,
//     // add other metadata here
//   }
// }
