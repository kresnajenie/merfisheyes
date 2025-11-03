import { redirect } from "next/navigation";

export default function SingleMoleculeViewerRedirect() {
  redirect("/?mode=sm");
}
