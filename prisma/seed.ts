import { PrismaClient } from "@prisma/client";

const prisma = new PrismaClient();

async function main() {
  const email = process.env.SUPER_ADMIN_EMAIL;

  if (!email) {
    console.log("SUPER_ADMIN_EMAIL not set — skipping super admin seed.");

    return;
  }

  const user = await prisma.user.findUnique({ where: { email } });

  if (!user) {
    console.log(
      `User with email "${email}" not found. Sign in first, then re-run the seed.`,
    );

    return;
  }

  if (user.role === "SUPER_ADMIN") {
    console.log(`User "${email}" is already SUPER_ADMIN.`);

    return;
  }

  await prisma.user.update({
    where: { email },
    data: { role: "SUPER_ADMIN" },
  });

  console.log(`Promoted "${email}" to SUPER_ADMIN.`);
}

main()
  .catch((e) => {
    console.error(e);
    process.exit(1);
  })
  .finally(() => prisma.$disconnect());
