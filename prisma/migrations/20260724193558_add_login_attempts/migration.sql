-- CreateTable
CREATE TABLE "login_attempts" (
    "id" TEXT NOT NULL,
    "identifier" TEXT NOT NULL,
    "kind" TEXT NOT NULL,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "login_attempts_pkey" PRIMARY KEY ("id")
);

-- CreateIndex
CREATE INDEX "login_attempts_identifier_kind_created_at_idx" ON "login_attempts"("identifier", "kind", "created_at");

