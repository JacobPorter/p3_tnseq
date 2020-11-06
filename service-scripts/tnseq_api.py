import os
import sys

LOG = sys.stderr

PatricUser = None


def authenticateByEnv(Session):
    if "KB_AUTH_TOKEN" in os.environ:
        LOG.write("reading auth key from environment\n")
        authenticateByString(os.environ.get('KB_AUTH_TOKEN'), Session)
        return True
    return authenticateByFile(None, Session)


def authenticateByString(tokenString, Session):
    Session.headers.update({'Authorization': tokenString})
    if "Authorization" in Session.headers:
        global PatricUser
        PatricUser = Session.headers["Authorization"].split(r"|")[3].split(
            "=")[1]
        LOG.write("Patric user = %s\n" % PatricUser)


def authenticateByFile(tokenFile=None, Session=None):
    if not tokenFile:
        tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
    if os.path.exists(tokenFile):
        LOG.write("reading auth key from file %s\n" % tokenFile)
        with open(tokenFile) as F:
            tokenString = F.read().rstrip()
            authenticateByString(tokenString, Session)
        return True
    return False
