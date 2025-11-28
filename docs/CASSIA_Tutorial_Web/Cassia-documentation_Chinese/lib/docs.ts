import fs from "fs"
import path from "path"
import matter from "gray-matter"

const docsDirectory = path.join(process.cwd(), "content/docs")

export function getDocBySlug(slug: string) {
  const realSlug = slug.replace(/\.md$/, "")
  const fullPath = path.join(docsDirectory, `${realSlug}.md`)

  try {
    const fileContents = fs.readFileSync(fullPath, "utf8")
    const { data, content } = matter(fileContents)

    // Fix the triple backtick escaping in markdown
    const processedContent = content
      // Replace escaped backticks with actual backticks
      .replace(/\\`\\`\\`/g, "```")
      // Make sure only content inside triple backticks is treated as code blocks
      .replace(/```(\s*)(r|R)?\s*\n([\s\S]*?)```/g, (match, space, lang, code) => {
        return `\`\`\`${lang || "r"}\n${code}\`\`\``
      })

    return {
      slug: realSlug,
      frontmatter: data,
      content: processedContent,
    }
  } catch (error) {
    console.error(`Error reading file ${fullPath}:`, error)
    return {
      slug: realSlug,
      frontmatter: { title: "Not Found" },
      content: "# Not Found\n\nThe requested document could not be found.",
    }
  }
}

export function getAllDocs() {
  try {
    const slugs = fs
      .readdirSync(docsDirectory)
      .filter((file) => file.endsWith(".md"))
      .map((file) => file.replace(/\.md$/, ""))

    const docs = slugs.map((slug) => getDocBySlug(slug))

    return docs
  } catch (error) {
    console.error("Error reading docs directory:", error)
    return []
  }
}
