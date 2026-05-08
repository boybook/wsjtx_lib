#include "wsjtx_lib.h"
#include "wsjtx_decode.h"
#include <cstring>
#include <string>
#include <fftw3.h>
#include <ctime>
#include <time.h>

static thread_local DataQueue<WsjtxMessage>* s_currentQueue = nullptr;
void wsjtx_set_message_queue(DataQueue<WsjtxMessage>* q) { s_currentQueue = q; }

extern "C" {

void wsjtx_decoded_(int *nutc, int *snr, float *dt, int *freq, char *decoded, int len)
{
	char message[38];
	std::strncpy(message, decoded, 37);
	message[37] = '\0';
	for (int i = 37; i != 0; i--) {
		if (message[i] == ' ' || message[i] == '\0') message[i] = '\0';
		else break;
	}
	if (!strstr(message, "DecodeFinished") && s_currentQueue) {
		s_currentQueue->push(WsjtxMessage(*nutc / 10000, (*nutc / 100) % 100, *nutc % 100, *snr, *dt, *freq, std::string(message)));
	}
}

void wsjtx_decoded_fst4_(int *nutc, float *sync, int *snr, float *dt, float *freq, char *decoded, int len)
{
	char message[38];
	std::strncpy(message, decoded, 37);
	message[37] = '\0';
	for (int i = 37; i != 0; i--) {
		if (message[i] == ' ' || message[i] == '\0') message[i] = '\0';
		else break;
	}
	if (!strstr(message, "DecodeFinished") && s_currentQueue) {
		s_currentQueue->push(WsjtxMessage(*nutc / 10000, (*nutc / 100) % 100, *nutc % 100, *snr, *sync, *dt, *freq, std::string(message)));
	}
}
}

void wstjx_decode::setDxInfo(const std::string& call, const std::string& grid) {
	dx_call_ = call;
	dx_grid_ = grid;
}
void wstjx_decode::setStationInfo(const std::string& myCall, const std::string& myGrid,
	const std::string& dxCall, const std::string& dxGrid) {
	my_call_ = myCall;
	my_grid_ = myGrid;
	dx_call_ = dxCall;
	dx_grid_ = dxGrid;
}
void wstjx_decode::setDecodeRange(int low, int high, int tol) { nfa_ = low; nfb_ = high; ntol_ = tol; }
void wstjx_decode::setDecodeControls(bool apDecode, int decodeDepth, int txFrequency, int qsoProgress) {
	ap_decode_ = apDecode;
	decode_depth_ = decodeDepth < 1 ? 1 : decodeDepth;
	tx_frequency_ = txFrequency;
	qso_progress_ = qsoProgress < 0 ? 0 : qsoProgress;
}

void wstjx_decode::decode(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int threads)
{
	samplebuffer.push(std::move(audiosamples));
	std::memset(&params, 0, sizeof(params));
	params.nmode = 8; params.ntrperiod = 60.0; params.nQSOProgress = qso_progress_;
	params.nfqso = freq; params.nftx = tx_frequency_;
	params.newdat = true; params.npts8 = 74736; params.nfa = nfa_;
	params.nfSplit = 2700; params.nfb = nfb_; params.ntol = ntol_;
	params.kin = 64800; params.nzhsym = 79; params.nsubmode = 0;
	params.nagain = false; params.ndepth = decode_depth_; params.lft8apon = ap_decode_;
	params.lapcqonly = false; params.ljt65apon = true; params.napwid = 75;
	params.ntxmode = 65; params.nmode = 8; params.minw = 0;
	params.nclearave = false; params.minSync = 0; params.emedelay = 0.0;
	params.dttol = 3; params.nlist = 0; params.listutc[0] = '\0';
	params.n2pass = 2; params.nranera = 6; params.naggressive = 0;
	params.nrobust = false; params.nexp_decode = 0;
	if (!my_call_.empty()) std::strncpy(params.mycall, my_call_.c_str(), sizeof(params.mycall) - 1);
	if (!my_grid_.empty()) std::strncpy(params.mygrid, my_grid_.c_str(), sizeof(params.mygrid) - 1);
	if (!dx_call_.empty()) std::strncpy(params.hiscall, dx_call_.c_str(), sizeof(params.hiscall) - 1);
	if (!dx_grid_.empty()) std::strncpy(params.hisgrid, dx_grid_.c_str(), sizeof(params.hisgrid) - 1);
	switch (mode) {
	case FT8: params.nmode = 8; break;
	case FT4: params.nmode = 5; break;
	default: return;
	}
	int nfsample = 12000;
	for (size_t i = 0; i < audiosamples.size(); i++)
		dec_data.d2[i] = (short int)(audiosamples[i] * 32768.0f);
	auto now = std::chrono::system_clock::now();
	time_t tt = std::chrono::system_clock::to_time_t(now);
	tm local_tm = *localtime(&tt);
	params.nutc = local_tm.tm_hour * 10000 + local_tm.tm_min * 100 + local_tm.tm_sec;
	fftwf_plan_with_nthreads(threads);
	multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
}

void wstjx_decode::decode(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int threads)
{
	std::memset(&params, 0, sizeof(params));
	params.nmode = 8; params.ntrperiod = 60.0; params.nQSOProgress = qso_progress_;
	params.nfqso = freq; params.nftx = tx_frequency_;
	params.newdat = true; params.npts8 = 74736; params.nfa = nfa_;
	params.nfSplit = 2700; params.nfb = nfb_; params.ntol = ntol_;
	params.kin = 64800; params.nzhsym = 50; params.nsubmode = 0;
	params.nagain = false; params.ndepth = decode_depth_; params.lft8apon = ap_decode_;
	params.lapcqonly = false; params.ljt65apon = true; params.napwid = 75;
	params.ntxmode = 65;
	switch (mode) {
	case FT8: params.nmode = 8; break;
	case FT4: params.nmode = 5; break;
	default: return;
	}
	params.minw = 0; params.nclearave = false; params.minSync = 0;
	params.emedelay = 0.0; params.dttol = 3; params.nlist = 0;
	params.listutc[0] = '\0'; params.n2pass = 2; params.nranera = 6;
	params.naggressive = 0; params.nrobust = false; params.nexp_decode = 0;
	if (!my_call_.empty()) std::strncpy(params.mycall, my_call_.c_str(), sizeof(params.mycall) - 1);
	if (!my_grid_.empty()) std::strncpy(params.mygrid, my_grid_.c_str(), sizeof(params.mygrid) - 1);
	if (!dx_call_.empty()) std::strncpy(params.hiscall, dx_call_.c_str(), sizeof(params.hiscall) - 1);
	if (!dx_grid_.empty()) std::strncpy(params.hisgrid, dx_grid_.c_str(), sizeof(params.hisgrid) - 1);
	int nfsample = 12000;
	for (size_t i = 0; i < audiosamples.size(); i++)
		dec_data.d2[i] = (short int)audiosamples[i];
	auto now = std::chrono::system_clock::now();
	time_t tt = std::chrono::system_clock::to_time_t(now);
	tm local_tm = *localtime(&tt);
	params.nutc = local_tm.tm_hour * 10000 + local_tm.tm_min * 100 + local_tm.tm_sec;
	fftwf_plan_with_nthreads(threads);
	multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
}
