import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import "./globals.css";
import { AuthProvider } from "@/components/auth/AuthProvider";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "CASSIA - Collective Agent System for Single-cell Interpretable Annotation",
  description: "AI-powered single-cell RNA-seq analysis platform for automated cell type annotation with comprehensive analysis reports. Built for biologists, powered by large language models.",
  keywords: ["single-cell", "RNA-seq", "cell annotation", "bioinformatics", "AI", "machine learning"],
  authors: [{ name: "CASSIA Team" }],
  openGraph: {
    title: "CASSIA - Collective Agent System for Single-cell Interpretable Annotation",
    description: "AI-powered single-cell RNA-seq analysis platform",
    type: "website",
  },
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" className="scroll-smooth">
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased min-h-screen bg-gradient-to-br from-slate-50 via-blue-50 to-indigo-100 dark:from-slate-900 dark:via-slate-800 dark:to-slate-900`}
      >
        <div className="relative min-h-screen">
          {/* Animated background pattern */}
          <div className="absolute inset-0 opacity-30">
            <div className="absolute top-0 left-0 w-96 h-96 bg-gradient-to-r from-blue-400 to-purple-500 rounded-full mix-blend-multiply filter blur-xl opacity-70 animate-float"></div>
            <div className="absolute top-0 right-0 w-96 h-96 bg-gradient-to-r from-purple-400 to-pink-500 rounded-full mix-blend-multiply filter blur-xl opacity-70 animate-float" style={{animationDelay: '2s'}}></div>
            <div className="absolute bottom-0 left-1/2 w-96 h-96 bg-gradient-to-r from-cyan-400 to-blue-500 rounded-full mix-blend-multiply filter blur-xl opacity-70 animate-float" style={{animationDelay: '4s'}}></div>
          </div>
          
          {/* Main content */}
          <div className="relative z-10">
            <AuthProvider>
              {children}
            </AuthProvider>
          </div>
        </div>
      </body>
    </html>
  );
}
