#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Frame Sync
# Author: Xiangyu Zhang
# Description: frame syncing
# GNU Radio version: 3.8.2.0

from distutils.version import StrictVersion

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print("Warning: failed to XInitThreads()")

from gnuradio import analog
from gnuradio import blocks
import pmt
from gnuradio import digital
from gnuradio import gr
from gnuradio.filter import firdes
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
import numpy

from gnuradio import qtgui

class frame_sync(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Frame Sync")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Frame Sync")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "frame_sync")

        try:
            if StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
                self.restoreGeometry(self.settings.value("geometry").toByteArray())
            else:
                self.restoreGeometry(self.settings.value("geometry"))
        except:
            pass

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 32000
        self.length_tag_key = length_tag_key = "# of bits"
        self.header_mod = header_mod = digital.constellation_bpsk()

        ##################################################
        # Blocks
        ##################################################
        self.digital_protocol_formatter_bb_0 = digital.protocol_formatter_bb(hdr_format, length_tag_key)
        self.digital_crc32_bb_0_0 = digital.crc32_bb(True, length_tag_key, True)
        self.digital_crc32_bb_0 = digital.crc32_bb(False, length_tag_key, True)
        self.digital_correlate_access_code_xx_ts_1 = digital.correlate_access_code_bb_ts(digital.packet_utils.default_access_code,
          1, length_tag_key)
        self.digital_constellation_decoder_cb_0 = digital.constellation_decoder_cb(header_mod.base())
        self.digital_chunks_to_symbols_xx_0_0_0_0_0_0_0_0 = digital.chunks_to_symbols_bc(header_mod.points(), 1)
        self.digital_chunks_to_symbols_xx_0_0_0_0_0 = digital.chunks_to_symbols_bc(header_mod.points(), 1)
        self.blocks_tagged_stream_to_pdu_0 = blocks.tagged_stream_to_pdu(blocks.byte_t, length_tag_key)
        self.blocks_tagged_stream_mux_0 = blocks.tagged_stream_mux(gr.sizeof_gr_complex*1, length_tag_key, 0)
        self.blocks_tag_gate_0_0 = blocks.tag_gate(gr.sizeof_gr_complex * 1, False)
        self.blocks_tag_gate_0_0.set_single_key("")
        self.blocks_repack_bits_bb_0_0_1_0 = blocks.repack_bits_bb(8, header_mod.bits_per_symbol(), length_tag_key, False, gr.GR_MSB_FIRST)
        self.blocks_repack_bits_bb_0_0_1 = blocks.repack_bits_bb(8, header_mod.bits_per_symbol(), length_tag_key, False, gr.GR_MSB_FIRST)
        self.blocks_repack_bits_bb_0_0_0 = blocks.repack_bits_bb(1, 8, length_tag_key, False, gr.GR_MSB_FIRST)
        self.blocks_pdu_to_tagged_stream_1 = blocks.pdu_to_tagged_stream(blocks.byte_t, length_tag_key)
        self.blocks_message_strobe_0_0_0_0 = blocks.message_strobe(pmt.cons(pmt.make_dict(), pmt.pmt_to_python.numpy_to_uvector(numpy.array([ord(c) for c in "Beyond 5G"], numpy.uint8))), 1000)
        self.blocks_message_debug_0 = blocks.message_debug()
        self.blocks_add_xx_0 = blocks.add_vcc(1)
        self.analog_noise_source_x_0 = analog.noise_source_c(analog.GR_GAUSSIAN, 0.7, 0)



        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.blocks_message_strobe_0_0_0_0, 'strobe'), (self.blocks_pdu_to_tagged_stream_1, 'pdus'))
        self.msg_connect((self.blocks_tagged_stream_to_pdu_0, 'pdus'), (self.blocks_message_debug_0, 'print'))
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.blocks_add_xx_0, 0), (self.blocks_tag_gate_0_0, 0))
        self.connect((self.blocks_pdu_to_tagged_stream_1, 0), (self.digital_crc32_bb_0, 0))
        self.connect((self.blocks_repack_bits_bb_0_0_0, 0), (self.digital_crc32_bb_0_0, 0))
        self.connect((self.blocks_repack_bits_bb_0_0_1, 0), (self.digital_chunks_to_symbols_xx_0_0_0_0_0, 0))
        self.connect((self.blocks_repack_bits_bb_0_0_1_0, 0), (self.digital_chunks_to_symbols_xx_0_0_0_0_0_0_0_0, 0))
        self.connect((self.blocks_tag_gate_0_0, 0), (self.digital_constellation_decoder_cb_0, 0))
        self.connect((self.blocks_tagged_stream_mux_0, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.digital_chunks_to_symbols_xx_0_0_0_0_0, 0), (self.blocks_tagged_stream_mux_0, 1))
        self.connect((self.digital_chunks_to_symbols_xx_0_0_0_0_0_0_0_0, 0), (self.blocks_tagged_stream_mux_0, 0))
        self.connect((self.digital_constellation_decoder_cb_0, 0), (self.digital_correlate_access_code_xx_ts_1, 0))
        self.connect((self.digital_correlate_access_code_xx_ts_1, 0), (self.blocks_repack_bits_bb_0_0_0, 0))
        self.connect((self.digital_crc32_bb_0, 0), (self.blocks_repack_bits_bb_0_0_1, 0))
        self.connect((self.digital_crc32_bb_0, 0), (self.digital_protocol_formatter_bb_0, 0))
        self.connect((self.digital_crc32_bb_0_0, 0), (self.blocks_tagged_stream_to_pdu_0, 0))
        self.connect((self.digital_protocol_formatter_bb_0, 0), (self.blocks_repack_bits_bb_0_0_1_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "frame_sync")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

    def get_length_tag_key(self):
        return self.length_tag_key

    def set_length_tag_key(self, length_tag_key):
        self.length_tag_key = length_tag_key

    def get_header_mod(self):
        return self.header_mod

    def set_header_mod(self, header_mod):
        self.header_mod = header_mod





def main(top_block_cls=frame_sync, options=None):

    if StrictVersion("4.5.0") <= StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    def quitting():
        tb.stop()
        tb.wait()

    qapp.aboutToQuit.connect(quitting)
    qapp.exec_()

if __name__ == '__main__':
    main()
