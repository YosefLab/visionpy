from flask import (
    Blueprint
)

bp = Blueprint('Signature', __name__, url_prefix='/Signature')

@bp.route('/Signature/Scores/<sig_name>', methods=('GET'))
def getScores(sig_name):
    '''
    '''
    if request.method == 'GET':
        return None

@bp.route('/Signature/Meta/<sig_name>', methods=('GET'))
def getMeta(sig_name):
    '''
    '''
    if request.method == 'GET':
        return None


@bp.route('/Signature/Info/<sig_name>', methods=('GET'))
def getInfo(sig_name):
    '''
    '''
    if request.method == 'GET':
        return None

@bp.route('/Signature/Expression/<sig_name>/<cluster_var>', methods=('GET'))
def getClusterExpression(sig_name, cluster_var):
    '''
    '''
    if request.method == 'GET':
        return None
