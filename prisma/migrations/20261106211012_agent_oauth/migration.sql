-- CreateTable
CREATE TABLE "oauth_clients" (
    "id" VARCHAR(1000) NOT NULL,
    "kind" VARCHAR(20) NOT NULL,
    "name" VARCHAR(200) NOT NULL,
    "uri" VARCHAR(1000),
    "logo_uri" VARCHAR(1000),
    "redirect_uris" TEXT[],
    "expires_at" TIMESTAMP(3),
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,
    "updated_at" TIMESTAMP(3) NOT NULL,

    CONSTRAINT "oauth_clients_pkey" PRIMARY KEY ("id")
);

-- CreateTable
CREATE TABLE "oauth_grants" (
    "id" TEXT NOT NULL,
    "user_id" TEXT NOT NULL,
    "client_id" VARCHAR(1000) NOT NULL,
    "scopes" TEXT[],
    "all_projects" BOOLEAN NOT NULL DEFAULT true,
    "project_ids" TEXT[] DEFAULT ARRAY[]::TEXT[],
    "resource" VARCHAR(1000) NOT NULL,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,
    "updated_at" TIMESTAMP(3) NOT NULL,
    "last_used_at" TIMESTAMP(3),
    "revoked_at" TIMESTAMP(3),

    CONSTRAINT "oauth_grants_pkey" PRIMARY KEY ("id")
);

-- CreateTable
CREATE TABLE "oauth_codes" (
    "code_hash" VARCHAR(64) NOT NULL,
    "grant_id" TEXT NOT NULL,
    "redirect_uri" VARCHAR(1000) NOT NULL,
    "code_challenge" VARCHAR(128) NOT NULL,
    "resource" VARCHAR(1000) NOT NULL,
    "scopes" TEXT[],
    "expires_at" TIMESTAMP(3) NOT NULL,
    "used_at" TIMESTAMP(3),
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "oauth_codes_pkey" PRIMARY KEY ("code_hash")
);

-- CreateTable
CREATE TABLE "oauth_tokens" (
    "id" TEXT NOT NULL,
    "grant_id" TEXT NOT NULL,
    "kind" VARCHAR(10) NOT NULL,
    "token_hash" VARCHAR(64) NOT NULL,
    "resource" VARCHAR(1000) NOT NULL,
    "scopes" TEXT[],
    "expires_at" TIMESTAMP(3) NOT NULL,
    "used_at" TIMESTAMP(3),
    "revoked_at" TIMESTAMP(3),
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "oauth_tokens_pkey" PRIMARY KEY ("id")
);

-- CreateIndex
CREATE UNIQUE INDEX "oauth_grants_user_id_client_id_key" ON "oauth_grants"("user_id", "client_id");

-- CreateIndex
CREATE INDEX "oauth_codes_grant_id_idx" ON "oauth_codes"("grant_id");

-- CreateIndex
CREATE UNIQUE INDEX "oauth_tokens_token_hash_key" ON "oauth_tokens"("token_hash");

-- CreateIndex
CREATE INDEX "oauth_tokens_grant_id_idx" ON "oauth_tokens"("grant_id");

-- AddForeignKey
ALTER TABLE "oauth_grants" ADD CONSTRAINT "oauth_grants_user_id_fkey" FOREIGN KEY ("user_id") REFERENCES "users"("id") ON DELETE CASCADE ON UPDATE CASCADE;

-- AddForeignKey
ALTER TABLE "oauth_grants" ADD CONSTRAINT "oauth_grants_client_id_fkey" FOREIGN KEY ("client_id") REFERENCES "oauth_clients"("id") ON DELETE CASCADE ON UPDATE CASCADE;

-- AddForeignKey
ALTER TABLE "oauth_codes" ADD CONSTRAINT "oauth_codes_grant_id_fkey" FOREIGN KEY ("grant_id") REFERENCES "oauth_grants"("id") ON DELETE CASCADE ON UPDATE CASCADE;

-- AddForeignKey
ALTER TABLE "oauth_tokens" ADD CONSTRAINT "oauth_tokens_grant_id_fkey" FOREIGN KEY ("grant_id") REFERENCES "oauth_grants"("id") ON DELETE CASCADE ON UPDATE CASCADE;
