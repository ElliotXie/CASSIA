import type { Metadata } from 'next';
import PipelineClient from './PipelineClient';

export const metadata: Metadata = {
  title: 'CASSIA Pipeline - Single-cell Analysis',
  description: 'Full CASSIA analysis pipeline for single-cell RNA-seq cell type annotation.',
};

export default function PipelinePage() {
  return <PipelineClient />;
}
