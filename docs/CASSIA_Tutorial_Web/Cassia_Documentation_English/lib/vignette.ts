import fs from "fs"
import path from "path"
import matter from "gray-matter"

const contentDirectory = path.join(process.cwd(), "content")

export type ProgrammingLang = "r" | "python"

export function getVignetteBySlug(slug: string, locale: string = 'en', programmingLang: ProgrammingLang = 'r') {
  const realSlug = slug.replace(/\.md$/, "")

  // Try locale-specific path with programming language first
  let fullPath = path.join(contentDirectory, locale, "vignette", programmingLang, `${realSlug}.md`)
  let isFallback = false
  let isPythonFallback = false

  // Fallback chain:
  // 1. Try requested locale + requested programming lang
  // 2. Try requested locale + R (if Python was requested but not found)
  // 3. Try English + requested programming lang
  // 4. Try English + R

  if (!fs.existsSync(fullPath)) {
    // If Python was requested but not found, try R version in same locale
    if (programmingLang === 'python') {
      const rPath = path.join(contentDirectory, locale, "vignette", "r", `${realSlug}.md`)
      if (fs.existsSync(rPath)) {
        fullPath = rPath
        isPythonFallback = true
      }
    }
  }

  if (!fs.existsSync(fullPath) && locale !== 'en') {
    // Try English version with requested programming lang
    const enPath = path.join(contentDirectory, 'en', "vignette", programmingLang, `${realSlug}.md`)
    if (fs.existsSync(enPath)) {
      fullPath = enPath
      isFallback = true
    } else if (programmingLang === 'python') {
      // Try English + R as final fallback
      const enRPath = path.join(contentDirectory, 'en', "vignette", "r", `${realSlug}.md`)
      if (fs.existsSync(enRPath)) {
        fullPath = enRPath
        isFallback = true
        isPythonFallback = true
      }
    }
  }

  try {
    const fileContents = fs.readFileSync(fullPath, "utf8")
    const { data, content } = matter(fileContents)

    // Fix the triple backtick escaping in markdown
    const processedContent = content
      // Replace escaped backticks with actual backticks
      .replace(/\\`\\`\\`/g, "```")

    return {
      slug: realSlug,
      frontmatter: data,
      content: processedContent,
      isFallback,
      isPythonFallback,
    }
  } catch (error) {
    console.error(`Error reading file ${fullPath}:`, error)
    return {
      slug: realSlug,
      frontmatter: { title: "Not Found" },
      content: "# Not Found\n\nThe requested vignette could not be found.",
      isFallback: false,
      isPythonFallback: false,
    }
  }
}

export function getAllVignettes(locale: string = 'en', programmingLang: ProgrammingLang = 'r') {
  // Get vignettes for the locale, with English fallback
  const vignetteDir = path.join(contentDirectory, 'en', 'vignette', 'r')

  try {
    if (!fs.existsSync(vignetteDir)) {
      return [];
    }

    // Use English R slugs as the canonical list
    const slugs = fs
      .readdirSync(vignetteDir)
      .filter((file) => file.endsWith(".md"))
      .map((file) => file.replace(/\.md$/, ""))

    const vignettes = slugs.map((slug) => getVignetteBySlug(slug, locale, programmingLang))

    return vignettes
  } catch (error) {
    console.error("Error reading vignette directory:", error)
    return []
  }
}

// Get all available programming languages for a given vignette
export function getAvailableProgrammingLanguages(slug: string, locale: string = 'en'): ProgrammingLang[] {
  const realSlug = slug.replace(/\.md$/, "")
  const available: ProgrammingLang[] = []

  const rPath = path.join(contentDirectory, locale, "vignette", "r", `${realSlug}.md`)
  const pythonPath = path.join(contentDirectory, locale, "vignette", "python", `${realSlug}.md`)

  if (fs.existsSync(rPath)) {
    available.push('r')
  }
  if (fs.existsSync(pythonPath)) {
    available.push('python')
  }

  // Fallback to English if nothing found
  if (available.length === 0 && locale !== 'en') {
    return getAvailableProgrammingLanguages(slug, 'en')
  }

  return available
}
