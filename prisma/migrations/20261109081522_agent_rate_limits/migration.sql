-- CreateTable
CREATE TABLE "rate_limit_counters" (
    "key" VARCHAR(300) NOT NULL,
    "window_start" TIMESTAMP(3) NOT NULL,
    "count" INTEGER NOT NULL DEFAULT 0,

    CONSTRAINT "rate_limit_counters_pkey" PRIMARY KEY ("key","window_start")
);

-- CreateIndex
CREATE INDEX "rate_limit_counters_window_start_idx" ON "rate_limit_counters"("window_start");
