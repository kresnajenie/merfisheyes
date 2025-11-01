import { redirect } from "next/navigation";

interface SingleMoleculeViewerIdPageProps {
  params: { id: string };
}

export default function SingleMoleculeViewerIdRedirect({
  params,
}: SingleMoleculeViewerIdPageProps) {
  const search = new URLSearchParams({ mode: "sm", dataset: params.id });

  redirect(`/?${search.toString()}`);
}
