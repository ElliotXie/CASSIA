import type { Metadata } from 'next';
import MergingClient from './MergingClient';

export const metadata: Metadata = {
  title: 'Annotation Merging Agent - CASSIA',
  description: 'Merge and group cell cluster annotations using AI to create broader cell type categories.',
};

export default function AnnotationMergingPage() {
  return <MergingClient />;
}
