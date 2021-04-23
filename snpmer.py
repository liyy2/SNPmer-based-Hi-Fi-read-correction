class SnpmerParser:
    def __init__(self, parsed_str):
        ret = parsed_str.split(':')
        self.snpmer_id = ret[0][:-1]
        self.reference = ret[0][-1]
        assert (self.reference == 'R') | (self. reference == 'A'), 'Expected reference to be R or A, but get {}'.format(self.reference)
        if ret[1][0] == '-':
            self.directionality = '-'
            ret[1] = ret[1:]
        else:
            self.directionality = '+'
        self.position_in_read = int(ret[1])
    
    def generate_identifier(self):
        return self.snpmer_id + self.reference



        