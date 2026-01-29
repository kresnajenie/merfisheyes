/**
 * Simple Database Connection Test
 *
 * Tests if DATABASE_URL and DIRECT_URL are configured correctly
 * and if the database schema exists.
 *
 * Usage:
 *   node scripts/test-db-connection.js
 */

const { execSync } = require('child_process');

// Colors
const colors = {
  reset: '\x1b[0m',
  green: '\x1b[32m',
  red: '\x1b[31m',
  yellow: '\x1b[33m',
  blue: '\x1b[34m',
  cyan: '\x1b[36m',
};

const log = {
  info: (msg) => console.log(`${colors.blue}ℹ${colors.reset} ${msg}`),
  success: (msg) => console.log(`${colors.green}✓${colors.reset} ${msg}`),
  error: (msg) => console.log(`${colors.red}✗${colors.reset} ${msg}`),
  warning: (msg) => console.log(`${colors.yellow}⚠${colors.reset} ${msg}`),
};

async function main() {
  console.log('\n' + colors.cyan + '🔍 Database Connection Test' + colors.reset + '\n');

  // Test 1: Check environment variables
  console.log(colors.cyan + '1. Checking Environment Variables' + colors.reset);

  const databaseUrl = process.env.DATABASE_URL;
  const directUrl = process.env.DIRECT_URL;

  if (!databaseUrl) {
    log.error('DATABASE_URL not set');
    console.log('');
    console.log('Set it in .env.local:');
    console.log('DATABASE_URL=postgresql://user:password@host:6543/postgres?pgbouncer=true');
    console.log('');
    process.exit(1);
  }
  log.success('DATABASE_URL is set');

  if (!directUrl) {
    log.warning('DIRECT_URL not set (optional for queries, required for migrations)');
  } else {
    log.success('DIRECT_URL is set');
  }

  // Mask password in output
  const maskedDbUrl = databaseUrl.replace(/:[^:@]+@/, ':***@');
  log.info(`Database: ${maskedDbUrl}`);
  console.log('');

  // Test 2: Check Prisma Client
  console.log(colors.cyan + '2. Checking Prisma Client' + colors.reset);

  try {
    execSync('npx prisma --version', { stdio: 'pipe' });
    log.success('Prisma CLI installed');
  } catch (error) {
    log.error('Prisma CLI not found');
    process.exit(1);
  }
  console.log('');

  // Test 3: Test database connection
  console.log(colors.cyan + '3. Testing Database Connection' + colors.reset);

  try {
    // Use the simpler approach - just try to validate the schema
    const output = execSync('npx prisma validate', {
      encoding: 'utf-8',
      stdio: 'pipe'
    });
    log.success('Prisma schema is valid');
  } catch (error) {
    log.error('Prisma validation failed');
    console.log('Error:', error.message);
  }

  // Try a simple connection test using Prisma
  try {
    const { PrismaClient } = require('@prisma/client');
    const prisma = new PrismaClient();

    // Simple connection test
    await prisma.$connect();
    log.success('Database connection successful');
    await prisma.$disconnect();
  } catch (error) {
    log.error('Database connection failed');
    console.log('');
    console.log('Error:', error.message);
    console.log('');
    console.log('Possible issues:');
    console.log('  - Wrong DATABASE_URL');
    console.log('  - Database not running');
    console.log('  - Network/firewall blocking connection');
    console.log('  - Wrong password');
    console.log('');
    process.exit(1);
  }
  console.log('');

  // Test 4: Check if schema exists
  console.log(colors.cyan + '4. Checking Database Schema' + colors.reset);

  try {
    const { PrismaClient } = require('@prisma/client');
    const prisma = new PrismaClient();

    // Query for tables
    const tables = await prisma.$queryRaw`
      SELECT tablename
      FROM pg_tables
      WHERE schemaname = 'public'
    `;

    await prisma.$disconnect();

    if (tables.length === 0) {
      log.warning('No tables found - database is empty');
      console.log('');
      console.log('Run migrations to create schema:');
      console.log('  npx prisma migrate deploy');
      console.log('');
    } else {
      log.success(`Found ${tables.length} table(s)`);

      // Check for expected tables
      const expectedTables = ['datasets', 'upload_sessions', 'upload_files', '_prisma_migrations'];
      const tableNames = tables.map(t => t.tablename);
      const foundExpected = expectedTables.filter(table => tableNames.includes(table));

      if (foundExpected.length === expectedTables.length) {
        log.success('All expected tables exist');
        console.log('');
        console.log('Tables found:');
        foundExpected.forEach(table => console.log(`  - ${table}`));
      } else {
        log.warning('Some expected tables are missing');
        console.log('');
        console.log('Expected:', expectedTables.join(', '));
        console.log('Found:', foundExpected.join(', ') || 'none');
        console.log('');
        console.log('Run migrations to create missing tables:');
        console.log('  npx prisma migrate deploy');
      }
    }
  } catch (error) {
    log.error('Failed to check schema');
    console.log('Error:', error.message);
  }
  console.log('');

  // Test 5: Check migration status
  console.log(colors.cyan + '5. Checking Migration Status' + colors.reset);

  try {
    const output = execSync('npx prisma migrate status', {
      encoding: 'utf-8',
      stdio: 'pipe'
    });

    if (output.includes('Database schema is up to date')) {
      log.success('All migrations applied');
    } else if (output.includes('following migration have not yet been applied')) {
      log.warning('Pending migrations detected');
      console.log('');
      console.log('Apply migrations:');
      console.log('  npx prisma migrate deploy');
    } else {
      log.info('Migration status:');
      console.log(output);
    }
  } catch (error) {
    // migrate status returns non-zero if migrations pending
    const output = error.stdout?.toString() || error.message;

    if (output.includes('following migration have not yet been applied')) {
      log.warning('Pending migrations detected');
      console.log('');
      console.log('Apply migrations:');
      console.log('  npx prisma migrate deploy');
    } else if (output.includes('Could not connect to the database')) {
      log.error('Cannot connect to database for migration check');
    } else {
      log.warning('Migration status check failed');
      console.log(output);
    }
  }
  console.log('');

  // Summary
  console.log(colors.cyan + '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━' + colors.reset);
  console.log(colors.green + '✓ Database connection test complete!' + colors.reset);
  console.log(colors.cyan + '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━' + colors.reset);
  console.log('');
}

main().catch(error => {
  console.error('Fatal error:', error);
  process.exit(1);
});
