#pragma once
#include <vector>
#include "DataBuffer.h"
#include "commons.h"
#include "fortran_interface.h"
#include "wsjtx_lib.h"

void wsjtx_set_message_queue(DataQueue<WsjtxMessage>* q);

class wstjx_decode
{
  public:
	void decode(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int threads = 1);
	void decode(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int threads = 1);
	void setDxInfo(const std::string& call, const std::string& grid);
	void setDecodeRange(int low, int high, int tol);
  private:
	params_t params;
	DataBuffer<float> samplebuffer;
	dec_data_t dec_data;
	int nfa_ = 200, nfb_ = 4000, ntol_ = 20;
	std::string dx_call_, dx_grid_;
};
