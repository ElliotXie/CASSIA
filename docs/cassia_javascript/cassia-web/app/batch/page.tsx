import type { Metadata } from 'next';
import BatchClient from './BatchClient';

export const metadata: Metadata = {
  title: 'Batch Processing - CASSIA',
  description: 'Parallel batch analysis for large single-cell RNA-seq datasets with CASSIA.',
};

export default function BatchPage() {
  return <BatchClient />;
}
