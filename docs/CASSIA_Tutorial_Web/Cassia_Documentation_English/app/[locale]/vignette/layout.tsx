import type React from "react"
import { Outline } from "@/components/outline"

export default function VignetteLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <div className="flex min-h-screen">
      <div className="flex-1 px-4 py-8 md:px-8 lg:px-12">
        <div className="mx-auto max-w-3xl">{children}</div>
      </div>
      <Outline />
    </div>
  )
} 