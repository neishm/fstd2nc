def run():
  from fstd2nc.__main__ import _fstd2nc_cmdline_trapped
  from cccbuffer import Buffer
  _fstd2nc_cmdline_trapped (buffer_type=Buffer)

if __name__ == '__main__':
  run()
