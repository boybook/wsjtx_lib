#include "wsjtx_lib.h"
#include "wsjtx_decode.h"
#include <cstring>
#include <string>
#include <fftw3.h>
#include <ctime>
#include <time.h>
#include <algorithm>

static thread_local DataQueue<WsjtxMessage>* s_currentQueue = nullptr;
static thread_local WsjtxDecodeStats s_last_stats{};
static thread_local int s_candidate_count = 0;
void wsjtx_set_message_queue(DataQueue<WsjtxMessage>* q) { s_currentQueue = q; }

extern "C" void wsjtx_decode_candidates_(int *candidates)
{
	if (candidates && *candidates > 0) s_candidate_count += *candidates;
}

extern "C" void wsjtx_decode_stats_(int *stage, int *candidates, int *decoded, int *average)
{
	s_last_stats.stage_symbols = stage ? *stage : 50;
	s_last_stats.candidate_count = (candidates && *candidates > 0) ? *candidates : s_candidate_count;
	s_last_stats.decoded_count = decoded ? *decoded : 0;
	s_last_stats.average_count = average ? *average : 0;
	s_candidate_count = 0;
}

WsjtxDecodeStats wsjtx_last_stats() { return s_last_stats; }

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

void wstjx_decode::setDecodeStage(const wsjtx_decode_stage_t& stage) {
	stage_symbols_ = stage.stage_symbols > 0 ? stage.stage_symbols : 50;
	slot_utc_ = stage.slot_utc;
	reset_state_ = stage.reset_state;
	nagain_ = stage.nagain;
	eme_delay_ms_ = std::max(stage.eme_delay_ms, 0);
	session_id_ = stage.session_id;
}

void wstjx_decode::decode(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int threads)
{
	s_candidate_count = 0;
	std::fill_n(dec_data.d2, 180000, static_cast<short int>(0));
	std::fill_n(dec_data.ss, 184 * NSMAX, 0.0f);
	const int stage = std::clamp(stage_symbols_, 41, 50);
	int nfsample = 12000;
	auto prepare = [&](int depth, int symbols, bool reset) {
		std::memset(&params, 0, sizeof(params));
		params.nmode = mode == FT4 ? 5 : 8;
		params.ntrperiod = mode == FT4 ? 7 : 15;
		params.nQSOProgress = qso_progress_;
		params.nfqso = freq; params.nftx = tx_frequency_;
		params.newdat = !nagain_; params.npts8 = symbols * 432; params.nfa = nfa_;
		params.nfSplit = 2700; params.nfb = nfb_; params.ntol = ntol_;
		params.kin = 64800; params.nzhsym = mode == FT8 ? symbols : 50; params.nsubmode = 0;
		params.nagain = nagain_; params.ndepth = depth; params.lft8apon = reset ? false : ap_decode_;
		params.lapcqonly = false; params.ljt65apon = true; params.napwid = 75;
		params.ntxmode = mode == FT4 ? 5 : 8; params.minw = 0;
		params.nclearave = false; params.minSync = 0; params.emedelay = eme_delay_ms_ / 1000.0f;
		params.dttol = 3; params.nlist = 0; params.listutc[0] = '\0';
		params.n2pass = 2; params.nranera = 6; params.naggressive = 0;
		params.nrobust = false; params.nexp_decode = 0;
		if (!my_call_.empty()) std::strncpy(params.mycall, my_call_.c_str(), sizeof(params.mycall) - 1);
		if (!my_grid_.empty()) std::strncpy(params.mygrid, my_grid_.c_str(), sizeof(params.mygrid) - 1);
		if (!dx_call_.empty()) std::strncpy(params.hiscall, dx_call_.c_str(), sizeof(params.hiscall) - 1);
		if (!dx_grid_.empty()) std::strncpy(params.hisgrid, dx_grid_.c_str(), sizeof(params.hisgrid) - 1);
		if (slot_utc_ >= 0) {
			params.nutc = slot_utc_;
		} else {
			auto now = std::chrono::system_clock::now();
			time_t tt = std::chrono::system_clock::to_time_t(now);
			tm local_tm = *localtime(&tt);
			params.nutc = local_tm.tm_hour * 10000 + local_tm.tm_min * 100 + local_tm.tm_sec;
		}
	};
	if (mode != FT8 && mode != FT4) return;
	fftwf_plan_with_nthreads(threads);
	auto loadAudio = [&]() {
		std::fill_n(dec_data.d2, 180000, static_cast<short int>(0));
		const size_t count = std::min<size_t>(audiosamples.size(), 180000);
		for (size_t i = 0; i < count; ++i) dec_data.d2[i] = static_cast<short int>(audiosamples[i] * 32768.0f);
	};
	if (reset_state_ && mode == FT8 && stage != 41) {
		loadAudio();
		prepare(decode_depth_, 41, true);
		multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
		if (stage > 47) {
			loadAudio();
			prepare(decode_depth_, 47, false);
			multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
		}
	}
	prepare(decode_depth_, mode == FT8 ? stage : 50, false);
	loadAudio();
	multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
}

void wstjx_decode::decode(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int threads)
{
	s_candidate_count = 0;
	std::fill_n(dec_data.d2, 180000, static_cast<short int>(0));
	std::fill_n(dec_data.ss, 184 * NSMAX, 0.0f);
	const int stage = std::clamp(stage_symbols_, 41, 50);
	int nfsample = 12000;
	auto prepare = [&](int depth, int symbols, bool reset) {
		std::memset(&params, 0, sizeof(params));
		params.nmode = mode == FT4 ? 5 : 8;
		params.ntrperiod = mode == FT4 ? 7 : 15;
		params.nQSOProgress = qso_progress_;
		params.nfqso = freq; params.nftx = tx_frequency_;
		params.newdat = !nagain_; params.npts8 = symbols * 432; params.nfa = nfa_;
		params.nfSplit = 2700; params.nfb = nfb_; params.ntol = ntol_;
		params.kin = 64800; params.nzhsym = mode == FT8 ? symbols : 50; params.nsubmode = 0;
		params.nagain = nagain_; params.ndepth = depth; params.lft8apon = reset ? false : ap_decode_;
		params.lapcqonly = false; params.ljt65apon = true; params.napwid = 75;
		params.ntxmode = mode == FT4 ? 5 : 8; params.minw = 0;
		params.nclearave = false; params.minSync = 0; params.emedelay = eme_delay_ms_ / 1000.0f;
		params.dttol = 3; params.nlist = 0; params.listutc[0] = '\0';
		params.n2pass = 2; params.nranera = 6; params.naggressive = 0;
		params.nrobust = false; params.nexp_decode = 0;
		if (!my_call_.empty()) std::strncpy(params.mycall, my_call_.c_str(), sizeof(params.mycall) - 1);
		if (!my_grid_.empty()) std::strncpy(params.mygrid, my_grid_.c_str(), sizeof(params.mygrid) - 1);
		if (!dx_call_.empty()) std::strncpy(params.hiscall, dx_call_.c_str(), sizeof(params.hiscall) - 1);
		if (!dx_grid_.empty()) std::strncpy(params.hisgrid, dx_grid_.c_str(), sizeof(params.hisgrid) - 1);
		if (slot_utc_ >= 0) params.nutc = slot_utc_;
		else {
			auto now = std::chrono::system_clock::now(); time_t tt = std::chrono::system_clock::to_time_t(now);
			tm local_tm = *localtime(&tt); params.nutc = local_tm.tm_hour * 10000 + local_tm.tm_min * 100 + local_tm.tm_sec;
		}
	};
	if (mode != FT8 && mode != FT4) return;
	fftwf_plan_with_nthreads(threads);
	auto loadAudio = [&]() {
		std::fill_n(dec_data.d2, 180000, static_cast<short int>(0));
		const size_t count = std::min<size_t>(audiosamples.size(), 180000);
		for (size_t i = 0; i < count; ++i) dec_data.d2[i] = audiosamples[i];
	};
	if (reset_state_ && mode == FT8 && stage != 41) {
		loadAudio();
		prepare(decode_depth_, 41, true);
		multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
		if (stage > 47) {
			loadAudio();
			prepare(decode_depth_, 47, false);
			multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
		}
	}
	prepare(decode_depth_, mode == FT8 ? stage : 50, false);
	loadAudio();
	multimode_decoder_(dec_data.ss, dec_data.d2, &params, &nfsample);
}
