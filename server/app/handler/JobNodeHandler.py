from database.Database import Database

class JobNodeHandler:

    @staticmethod
    async def Register(node_uniq, name, priority):
        node = await Database.JobNodeLoad(node_uniq)
        if node is None:
            node = await Database.JobNodeCreate(node_uniq, name, priority)
        else:
            await node.Refresh(name, priority)
        return node

    @staticmethod
    async def GetAll():
        return await Database.JobNodeLoadAll()

    @staticmethod
    async def Get(node_uniq):
        node = await Database.JobNodeLoad(node_uniq)
        if node is None:
            raise Exception('JobNode not found')
        return node
