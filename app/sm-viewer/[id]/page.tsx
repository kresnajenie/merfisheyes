import { redirect } from "next/navigation";

interface SingleMoleculeViewerIdPageProps {
  params: Promise<{ id: string }>;
}

export default async function SingleMoleculeViewerIdRedirect({
  params,
}: SingleMoleculeViewerIdPageProps) {
  const { id } = await params;
  const search = new URLSearchParams({ mode: "sm", dataset: id });

  redirect(`/?${search.toString()}`);
}
