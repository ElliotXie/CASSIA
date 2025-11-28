import { getVignetteBySlug, getAllVignettes } from "@/lib/vignette"
import DocPageClient from "../../docs/[slug]/DocPageClient"
import { notFound } from "next/navigation"

interface VignettePageProps {
  params: {
    slug: string
  }
}

// This is the Server Component for Vignette pages
export default async function VignettePage({ params }: VignettePageProps) {
  // Fetch data on the server
  const vignette = getVignetteBySlug(params.slug)

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
    content: vignette.content
  }

  // Pass the fetched data and params to the DocPageClient
  return <DocPageClient doc={processedVignette} params={params} />
}

// Add this function to generate static paths at build time
export async function generateStaticParams() {
  const vignettes = getAllVignettes()
  return vignettes.map((vignette) => ({
    slug: vignette.slug,
  }))
}

// Optional: Add generateMetadata to set page title, etc.
// import { Metadata } from 'next'
// export async function generateMetadata({ params }: VignettePageProps): Promise<Metadata> {
//   const vignette = getVignetteBySlug(params.slug)
//   if (!vignette) {
//     return { title: "Not Found" }
//   }
//   return {
//     title: vignette.frontmatter.title,
//     description: vignette.frontmatter.description || `Vignette for ${vignette.frontmatter.title}`,
//   }
// } 