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
  title: "Full Workflow Best Practices for Single-Cell RNA-seq Analysis",
  description: "Comprehensive guide covering quality control, doublet removal, and background determination for scRNA-seq analysis."
}

// This is the Server Component for the Full Workflow Best Practices page
export default function FullWorkflowBestPracticesPage() {
  // Hardcode the slug since this is a specific page
  const slug = "full-workflow-best-practices"
  
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
      title: vignette.frontmatter.title || "Full Workflow Best Practices for Single-Cell RNA-seq Analysis"
    },
    content: vignette.content
  }

  // Pass the fetched data to the DocPageClient
  return <DocPageClient doc={processedVignette} params={{ slug }} />
} 