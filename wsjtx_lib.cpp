#include "wsjtx_lib.h"
#include "wsjtx_decode.h"
#include "wsjtx_encode.h"
#include "constants.h"
#include <memory>
#include <fftw3.h>

static int s_Test = 0;
int wsjtx_libTest() { return ++s_Test; }

wsjtx_lib::wsjtx_lib() { fftwf_init_threads(); }

void wsjtx_lib::setDxCall(const std::string& call) { dx_call_ = call; }
void wsjtx_lib::setDxGrid(const std::string& grid) { dx_grid_ = grid; }

bool wsjtx_lib::pullMessage(WsjtxMessage &msg) { return messageQueue_.pull(msg); }

void wsjtx_lib::decode(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int thread)
{
	std::unique_ptr<wstjx_decode> ptr = std::make_unique<wstjx_decode>();
	if (!dx_call_.empty() || !dx_grid_.empty()) ptr->setDxInfo(dx_call_, dx_grid_);
	wsjtx_set_message_queue(&messageQueue_);
	ptr->decode(mode, audiosamples, freq, thread);
	wsjtx_set_message_queue(nullptr);
}

void wsjtx_lib::decode(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int thread)
{
	std::unique_ptr<wstjx_decode> ptr = std::make_unique<wstjx_decode>();
	if (!dx_call_.empty() || !dx_grid_.empty()) ptr->setDxInfo(dx_call_, dx_grid_);
	wsjtx_set_message_queue(&messageQueue_);
	ptr->decode(mode, audiosamples, freq, thread);
	wsjtx_set_message_queue(nullptr);
}

int wspr_decode(std::vector<std::complex<float>> &iqdat, int samples,
	decoder_options options, std::vector<struct decoder_results> &decodes, int threads);

std::vector<struct decoder_results> wsjtx_lib::wspr_decode(WsjtxIQSampleVector &iqsignal, decoder_options options)
{
	std::vector<decoder_results> results;
	::wspr_decode(iqsignal, iqsignal.size(), options, results, 4);
	return results;
}

std::vector<float> wsjtx_lib::encode(wsjtxMode mode, int frequency, std::string message, std::string &messagesend)
{
	switch (mode) {
	case FT8: {
		auto ptr = std::make_unique<wsjtx_encode>();
		return ptr->encode_ft8(mode, frequency, message, messagesend);
	}
	case FT4: {
		auto ptr = std::make_unique<wsjtx_encode>();
		return ptr->encode_ft4(mode, frequency, message, messagesend);
	}
	default: return {};
	}
}
