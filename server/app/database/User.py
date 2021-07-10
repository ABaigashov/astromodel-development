import time

class User:
    pool = None

    @staticmethod
    async def Init(pool):
        User.pool = pool
        async with pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    CREATE SEQUENCE IF NOT EXISTS seq_users_uid
                ''')
                await cur.execute('''
                    CREATE TABLE IF NOT EXISTS users (
                        uid BIGINT PRIMARY KEY DEFAULT nextval('seq_users_uid'),
                        external_uid BIGINT NOT NULL,
                        login VARCHAR(32) NOT NULL,
                        name VARCHAR(64),
                        session VARCHAR(64),
                        date_created BIGINT,
                        date_expired BIGINT
                    )
                ''')
                await cur.execute('''
                    CREATE INDEX IF NOT EXISTS idx_users_uid
                    ON users (uid)
                ''')
                await cur.execute('''
                    CREATE INDEX IF NOT EXISTS idx_users_external_uid
                    ON users (external_uid)
                ''')
                await cur.execute('''
                    CREATE INDEX IF NOT EXISTS idx_users_login
                    ON users (login)
                ''')
                await cur.execute('''
                    CREATE INDEX IF NOT EXISTS idx_users_session
                    ON users (session)
                ''')
                await cur.execute('''
                    CREATE INDEX IF NOT EXISTS idx_users_created
                    ON users (date_created)
                ''')
                await cur.execute('''
                    CREATE INDEX IF NOT EXISTS idx_users_expired
                    ON users (date_expired)
                ''')

    @staticmethod
    async def Create(external_uid, login, name):
        user = User()
        user.external_uid = external_uid
        user.login = login
        user.name = name
        user.session = ''
        user.date_created = int(time.time())
        user.date_expired = 0

        async with User.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    INSERT INTO users (
                        external_uid,
                        login,
                        name,
                        session,
                        date_created,
                        date_expired
                    )
                    VALUES ({}, '{}','{}','{}',{},{})
                    RETURNING uid
                '''.format
                    (
                    user.external_uid,
                    user.login,
                    user.name,
                    user.session,
                    user.date_created,
                    user.date_expired
                    )
                )

                result = await cur.fetchone()
                user.uid = result[0]

        return user

    @staticmethod
    async def Update(user):
        async with User.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    UPDATE users SET
                        name='{}',
                        session='{}',
                        date_created={},
                        date_expired={}
                    WHERE uid={}
                '''.format
                    (
                    user.name,
                    user.session,
                    user.date_created,
                    user.date_expired,
                    user.uid
                    )
                )
        return user

    @staticmethod
    async def Load(uid):
        async with User.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    SELECT * FROM users WHERE uid='{}' LIMIT 1'''.format(uid)
                )
                try:
                    result = await cur.fetchone()
                    user = User()
                    user.uid = result[0]
                    user.external_uid = result[1]
                    user.login = result[2]
                    user.name = result[3]
                    user.session = result[4]
                    user.date_created = result[5]
                    user.date_expired = result[6]
                    return user
                except Exception as e:
                    return None

    @staticmethod
    async def LoadByExternalUid(external_uid):
        async with User.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    SELECT * FROM users WHERE external_uid={} LIMIT 1'''.format(external_uid)
                )
                try:
                    result = await cur.fetchone()
                    user = User()
                    user.uid = result[0]
                    user.external_uid = result[1]
                    user.login = result[2]
                    user.name = result[3]
                    user.session = result[4]
                    user.date_created = result[5]
                    user.date_expired = result[6]
                    return user
                except Exception as e:
                    return None

    @staticmethod
    async def LoadByLogin(login):
        async with User.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    SELECT * FROM users WHERE login='{}' LIMIT 1'''.format(login)
                )
                try:
                    result = await cur.fetchone()
                    user = User()
                    user.uid = result[0]
                    user.external_uid = result[1]
                    user.login = result[2]
                    user.name = result[3]
                    user.session = result[4]
                    user.date_created = result[5]
                    user.date_expired = result[6]
                    return user
                except Exception as e:
                    return None

    @staticmethod
    async def LoadBySession(session):
        async with User.pool.acquire() as conn:
            async with conn.cursor() as cur:
                await cur.execute('''
                    SELECT * FROM users WHERE session='{}' LIMIT 1'''.format(session)
                )
                try:
                    result = await cur.fetchone()
                    user = User()
                    user.uid = result[0]
                    user.external_uid = result[1]
                    user.login = result[2]
                    user.name = result[3]
                    user.session = result[4]
                    user.date_created = result[5]
                    user.date_expired = result[6]
                    return user
                except Exception as e:
                    return None

    async def Authorize(self, session, duration):
        self.session = session
        self.date_expired = int(time.time()) + int(duration)
        await User.Update(self)

    async def Refresh(self, name):
        self.name = name
        await User.Update(self)

    def IsAuthorized(self):
        return self.date_expired > int(time.time())
