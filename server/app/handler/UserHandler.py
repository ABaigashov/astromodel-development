from database.Database import Database

class UserHandler:

    @staticmethod
    async def Register(external_uid, login, name):
        user = await Database.UserLoadByExternalUid(external_uid)
        if user is None:
            user = await Database.UserCreate(external_uid, login, name)
        else:
            await user.Refresh(name)

        return user

    async def Refresh(session, duration):
        user = await Database.UserLoadBySession(session)
        if user is None:
            raise Exception('User not found')
        if not user.IsAuthorized():
            raise Exception('User not authorized')

        await user.Authorize(session, duration);
        return user

    @staticmethod
    async def Authorize(external_uid, session, duration):
        user = await Database.UserLoadByExternalUid(external_uid)
        if user is None:
            raise Exception('User not found')

        await user.Authorize(session, duration);
        return user
