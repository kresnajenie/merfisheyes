-- Whose fault a FAILED processing run was. Auto-set from the error text on
-- failure; admin-overridable from the processing dashboard.
CREATE TYPE "DatasetFaultCategory" AS ENUM ('USER', 'PLATFORM', 'UNKNOWN');

ALTER TABLE "datasets" ADD COLUMN "fault_category" "DatasetFaultCategory";
