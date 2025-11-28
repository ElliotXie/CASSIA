import { getVignetteBySlug } from "@/lib/vignette"
import DocPageClient from "../../docs/[slug]/DocPageClient"
import { notFound } from "next/navigation"
import { Metadata } from "next"

interface VignettePageProps {
  params: {
    slug: string
  }
}

export const metadata: Metadata = {
  title: "Clustering and Annotation with Seurat",
  description: "Learn how to perform clustering and annotation starting from a raw Seurat object."
}

// This is the Server Component for the Clustering and Annotation page
export default function ClusteringAndAnnotationPage() {
  // Hardcode the slug since this is a specific page
  const slug = "clustering-and-annotation"
  
  // Fetch data on the server
  const vignette = getVignetteBySlug(slug)

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
      title: vignette.frontmatter.title || "Clustering and Annotation with Seurat"
    },
    content: vignette.content
  }

  // Pass the fetched data to the DocPageClient
  return <DocPageClient doc={processedVignette} params={{ slug }} />
} 