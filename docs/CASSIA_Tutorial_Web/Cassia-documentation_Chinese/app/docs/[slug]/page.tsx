import { getDocBySlug } from "@/lib/docs"
import DocPageClient from "./DocPageClient"
import { notFound } from "next/navigation"

interface DocPageProps {
  params: {
    slug: string
  }
}

// This is the Server Component
export default function DocPage({ params }: DocPageProps) {
  // Fetch data on the server
  const doc = getDocBySlug(params.slug)

  // Handle case where doc is not found
  if (!doc) {
    notFound()
  }

  // Pass the fetched data and params to the Client Component
  return <DocPageClient doc={doc} params={params} />
}

// Optional: Add generateStaticParams if you want to pre-render these pages at build time
// import { getAllDocs } from "@/lib/docs"
// export async function generateStaticParams() {
//  const docs = getAllDocs(["slug"])
//  return docs.map((doc) => ({
//    slug: doc.slug,
//  }))
// }

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
