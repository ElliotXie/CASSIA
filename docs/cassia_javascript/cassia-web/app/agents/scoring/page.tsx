import type { Metadata } from 'next';
import ScoringClient from './ScoringClient';

export const metadata: Metadata = {
  title: 'Scoring Agent - CASSIA',
  description: 'Score and evaluate cell type annotations using AI-powered analysis.',
};

export default function ScoringAgentPage() {
  return <ScoringClient />;
}
