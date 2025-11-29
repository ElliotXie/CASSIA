import fs from "fs"
import path from "path"
import matter from "gray-matter"

const contentDirectory = path.join(process.cwd(), "content")

export type ProgrammingLang = "r" | "python"

export function getDocBySlug(slug: string, locale: string = 'en', programmingLang: ProgrammingLang = 'r') {
  const realSlug = slug.replace(/\.md$/, "")

  // Try locale-specific path with programming language first
  let fullPath = path.join(contentDirectory, locale, "docs", programmingLang, `${realSlug}.md`)
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
      const rPath = path.join(contentDirectory, locale, "docs", "r", `${realSlug}.md`)
      if (fs.existsSync(rPath)) {
        fullPath = rPath
        isPythonFallback = true
      }
    }
  }

  if (!fs.existsSync(fullPath) && locale !== 'en') {
    // Try English version with requested programming lang
    const enPath = path.join(contentDirectory, 'en', "docs", programmingLang, `${realSlug}.md`)
    if (fs.existsSync(enPath)) {
      fullPath = enPath
      isFallback = true
    } else if (programmingLang === 'python') {
      // Try English + R as final fallback
      const enRPath = path.join(contentDirectory, 'en', "docs", "r", `${realSlug}.md`)
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
      // Make sure only content inside triple backticks is treated as code blocks
      .replace(/```(\s*)(r|R|python|Python)?\s*\n([\s\S]*?)```/g, (match, space, lang, code) => {
        return `\`\`\`${lang || "r"}\n${code}\`\`\``
      })

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
      content: "# Not Found\n\nThe requested document could not be found.",
      isFallback: false,
      isPythonFallback: false,
    }
  }
}

export function getAllDocs(locale: string = 'en', programmingLang: ProgrammingLang = 'r') {
  // Get docs for the locale, with English fallback
  const docsDir = path.join(contentDirectory, 'en', 'docs', 'r')

  try {
    // Use English R slugs as the canonical list
    const slugs = fs
      .readdirSync(docsDir)
      .filter((file) => file.endsWith(".md"))
      .map((file) => file.replace(/\.md$/, ""))

    const docs = slugs.map((slug) => getDocBySlug(slug, locale, programmingLang))

    return docs
  } catch (error) {
    console.error("Error reading docs directory:", error)
    return []
  }
}

// Get all available programming languages for a given doc
export function getAvailableProgrammingLanguages(slug: string, locale: string = 'en'): ProgrammingLang[] {
  const realSlug = slug.replace(/\.md$/, "")
  const available: ProgrammingLang[] = []

  const rPath = path.join(contentDirectory, locale, "docs", "r", `${realSlug}.md`)
  const pythonPath = path.join(contentDirectory, locale, "docs", "python", `${realSlug}.md`)

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
