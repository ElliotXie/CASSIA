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
  title: "Basic Annotation with Marker Files",
  description: "Learn how to use CASSIA for cell type annotation when you already have marker gene lists prepared."
}

// This is the Server Component for the Basic Annotation page
export default function BasicAnnotationPage() {
  // Hardcode the slug since this is a specific page
  const slug = "basic-annotation"
  
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
      title: vignette.frontmatter.title || "Basic Annotation with Marker Files"
    },
    content: vignette.content
  }

  // Pass the fetched data to the DocPageClient
  return <DocPageClient doc={processedVignette} params={{ slug }} />
} 