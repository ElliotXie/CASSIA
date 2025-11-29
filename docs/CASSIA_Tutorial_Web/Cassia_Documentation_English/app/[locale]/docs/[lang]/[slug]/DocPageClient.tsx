"use client"

import ReactMarkdown, { Components } from "react-markdown"
import rehypeSlug from "rehype-slug"
import rehypeAutolinkHeadings from "rehype-autolink-headings"
import { useState } from "react"
import { notFound } from "next/navigation"
import { ClassAttributes, HTMLAttributes } from "react"
import { useTranslations } from "next-intl"
import { type ProgrammingLang } from "@/lib/docs"

// Type for the doc object
type DocType = {
  slug: string
  frontmatter: { [key: string]: any; title: string }
  content: string
  isFallback?: boolean
  isPythonFallback?: boolean
}

// Helper function to create slug IDs for headings
function slugify(text: string) {
  return text
    .toString()
    .toLowerCase()
    .replace(/\s+/g, "-")
    .replace(/[^\w-]+/g, "")
    .replace(/--+/g, "-")
    .replace(/^-+/, "")
    .replace(/-+$/, "")
}

interface DocPageClientProps {
  params: {
    slug: string
    locale: string
    lang: ProgrammingLang
  }
  doc: DocType
}

// Define type for code component props provided by ReactMarkdown
interface CodeProps extends ClassAttributes<HTMLElement>, HTMLAttributes<HTMLElement> {
  node?: any
  inline?: boolean
  className?: string
  children?: React.ReactNode
}

export default function DocPageClient({ params, doc }: DocPageClientProps) {
  const [copiedId, setCopiedId] = useState<string | null>(null)
  const t = useTranslations("common")
  const tProg = useTranslations("programmingLanguage")

  if (!doc) {
    console.error("Doc not found in client component, should have been caught by server component.")
    notFound()
  }

  const handleCopy = (text: string, id: string) => {
    navigator.clipboard
      .writeText(text)
      .then(() => {
        setCopiedId(id)
        setTimeout(() => setCopiedId(null), 2000)
      })
      .catch((err) => {
        console.error("Failed to copy: ", err)
      })
  }

  // Define components object with proper typing for code
  const markdownComponents: Components = {
    p: ({ node, children }) => {
      if (node && node.children[0] && node.children[0].type === 'element' && node.children[0].tagName === 'div') {
        return <>{children}</>
      }
      return <p>{children}</p>
    },
    code: ({ node, inline, className, children, ...props }: CodeProps) => {
      const isInline = !node?.position?.start.line === node?.position?.end.line || !className?.startsWith("language-");

      if (isInline) {
        return (
          <code className="bg-muted px-1.5 py-0.5 rounded text-sm font-medium text-foreground" {...props}>
            {children}
          </code>
        )
      }

      const match = /language-(\w+)/.exec(className || "")
      const language = match ? match[1] : ""
      const content = String(children)
      const codeId = `code-${content.length}-${content.slice(0, 20).replace(/\W/g, "")}`
      const isCopied = copiedId === codeId

      // Format language display
      const displayLang = language === "r" || language === "R" ? "R" :
                          language === "python" || language === "Python" ? "Python" : language

      return (
        <div className="relative my-4">
          <pre className="rounded-md bg-slate-50 dark:bg-slate-800 p-4 overflow-x-auto">
            <code className={`language-${language || "text"} text-foreground font-medium`} {...props}>
              {children}
            </code>
          </pre>
          {language && (
            <div className="absolute top-2 left-3 text-xs font-medium text-slate-500 dark:text-slate-400">
              {displayLang}
            </div>
          )}
          <button
            onClick={() => handleCopy(content, codeId)}
            className="absolute top-2 right-2 rounded bg-slate-200/80 dark:bg-slate-700/80 p-1 text-xs font-medium text-slate-600 hover:bg-slate-300 dark:text-slate-300 dark:hover:bg-slate-600 transition-colors"
            aria-label={isCopied ? "Copied" : "Copy code"}
          >
            {isCopied ? (
              <svg
                xmlns="http://www.w3.org/2000/svg"
                width="14"
                height="14"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              >
                <path d="M20 6L9 17l-5-5"></path>
              </svg>
            ) : (
              <svg
                xmlns="http://www.w3.org/2000/svg"
                width="14"
                height="14"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              >
                <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
                <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
              </svg>
            )}
          </button>
        </div>
      )
    },
    h1: ({ node, ...props }) => <h1 id={slugify(props.children as string)} {...props} />,
    h2: ({ node, ...props }) => <h2 id={slugify(props.children as string)} {...props} />,
    h3: ({ node, ...props }) => <h3 id={slugify(props.children as string)} {...props} />,
    h4: ({ node, ...props }) => <h4 id={slugify(props.children as string)} {...props} />,
    h5: ({ node, ...props }) => <h5 id={slugify(props.children as string)} {...props} />,
    h6: ({ node, ...props }) => <h6 id={slugify(props.children as string)} {...props} />,
    img: ({ node, ...props }) => {
      const altTextParts = (props.alt || "").split('|');
      const altText = altTextParts[0];

      let imgStyle: { [key: string]: string } = {};

      if (altTextParts.length > 1) {
        const sizeInfo = altTextParts[1];
        const sizeParts = sizeInfo.split(',');

        sizeParts.forEach(part => {
          const [key, value] = part.split('=');
          if (key && value) {
            imgStyle[key.trim()] = value.trim();
          }
        });
      }

      return (
        <img
          {...props}
          alt={altText}
          style={imgStyle}
          className="rounded-lg my-6 max-w-full"
          loading="lazy"
        />
      );
    },
  }

  return (
    <article className="prose prose-slate dark:prose-invert max-w-none">
      {/* Show fallback notice for locale */}
      {doc.isFallback && (
        <div className="mb-4 p-3 bg-yellow-100 dark:bg-yellow-900/30 border border-yellow-300 dark:border-yellow-700 rounded-md text-sm text-yellow-800 dark:text-yellow-200">
          {t("translationNotAvailable")}
        </div>
      )}
      {/* Show "Coming Soon" notice for Python fallback */}
      {doc.isPythonFallback && (
        <div className="mb-4 p-3 bg-blue-100 dark:bg-blue-900/30 border border-blue-300 dark:border-blue-700 rounded-md text-sm text-blue-800 dark:text-blue-200">
          {tProg("pythonComingSoon")}
        </div>
      )}
      <h1>{doc.frontmatter.title}</h1>
      <ReactMarkdown
        rehypePlugins={[
          rehypeSlug,
          [
            rehypeAutolinkHeadings,
            {
              behavior: "wrap",
              properties: {
                className: ["anchor-link"],
              },
            },
          ],
        ]}
        components={markdownComponents}
      >
        {doc.content}
      </ReactMarkdown>
    </article>
  )
}
