class Logger(object):
    def __init__(self):
        self._logger = None
    
    def init(self, logdir, name = 'log'):
        if self._logger is None:
            import logging
            if not os.path.exists(logdir):
                os.makedirs(logdir)
            log_file = os.path.join(logdir, name)
            if os.path.exists(log_file):
                os.remove(log_file)
            self._logger = logging.getLogger()
            self._logger.setLevel('INFO')
            fh = logging.FileHandler(log_file)
            ch = logging.StreamHandler()
            self._logger.addHandler(fh)
            self._logger.addHandler(ch)
        
    def info(self, str_info):
        now = datetime.now()
        display_now = str(now).split(' ')[1][:-3]
        self.init('../log', 'tem.log')
        self._logger.info('[' + display_now + ']' + ' ' + str_info)

mylogger = Logger()    
printf = mylogger.info



    