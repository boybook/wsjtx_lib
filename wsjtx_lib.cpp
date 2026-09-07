#include "wsjtx_lib.h"
#include "wsjtx_decode.h"
#include "wsjtx_encode.h"
#include "constants.h"
#include <memory>
#include <fftw3.h>
#include <stdexcept>
#include <algorithm>

static int s_Test = 0;
int wsjtx_libTest() { return ++s_Test; }

wsjtx_lib::wsjtx_lib() { fftwf_init_threads(); }

void wsjtx_lib::setDxCall(const std::string& call) { std::scoped_lock lock(decodeMutex_); dx_call_ = call; }
void wsjtx_lib::setDxGrid(const std::string& grid) { std::scoped_lock lock(decodeMutex_); dx_grid_ = grid; }

void wsjtx_lib::setDecodeRange(int lowFreq, int highFreq, int tolerance)
{
	std::scoped_lock lock(decodeMutex_);
	decode_low_  = lowFreq;
	decode_high_ = highFreq;
	decode_tol_  = tolerance;
}

void wsjtx_lib::setDecodeStationInfo(const std::string& myCall, const std::string& myGrid,
	const std::string& dxCall, const std::string& dxGrid)
{
	std::scoped_lock lock(decodeMutex_);
	my_call_ = myCall;
	my_grid_ = myGrid;
	dx_call_ = dxCall;
	dx_grid_ = dxGrid;
}

void wsjtx_lib::setDecodeControls(bool apDecode, int decodeDepth, int txFrequency, int qsoProgress)
{
	std::scoped_lock lock(decodeMutex_);
	ap_decode_ = apDecode;
	decode_depth_ = decodeDepth < 1 ? 1 : decodeDepth;
	tx_frequency_ = txFrequency;
	qso_progress_ = qsoProgress < 0 ? 0 : qsoProgress;
}

void wsjtx_lib::setDecodeStage(const wsjtx_decode_stage_t& stage)
{
	std::scoped_lock lock(decodeMutex_);
	setDecodeStageLocked(stage);
}

void wsjtx_lib::setDecodeStageLocked(const wsjtx_decode_stage_t& stage)
{
	if (stage.stage_symbols > 0 && stage.stage_symbols != 41 && stage.stage_symbols != 47
		&& stage.stage_symbols != 49 && stage.stage_symbols != 50) {
		throw std::invalid_argument("FT8 decode stage must be 41, 47, 49, or 50");
	}
	const bool newSession = !stage.session_id.empty() && stage.session_id != activeSessionId_;
	if (newSession) {
		// A new slot owns a fresh result set. Discard messages left by a legacy
		// caller before resetting the native decoder lifecycle.
		WsjtxMessage stale;
		while (messageQueue_.pull(stale)) {}
		activeSessionId_ = stage.session_id;
		activeSessionStage_ = 0;
		activeSessionUtc_ = stage.slot_utc;
		resetState_ = true;
		sessionConfigured_ = true;
	}
	if (!stage.session_id.empty() && stage.stage_symbols > 0 && stage.stage_symbols < activeSessionStage_) {
		throw std::invalid_argument("decode stages must be monotonic within a session");
	}
	stageSymbols_ = stage.stage_symbols > 0 ? stage.stage_symbols : 50;
	slotUtc_ = stage.slot_utc;
	nagain_ = stage.nagain;
	emeDelayMs_ = std::max(stage.eme_delay_ms, 0);
	if (!stage.session_id.empty()) {
		activeSessionStage_ = std::max(activeSessionStage_, stage.stage_symbols);
		activeSessionUtc_ = stage.slot_utc;
	} else {
		activeSessionId_.clear();
		activeSessionStage_ = 0;
		activeSessionUtc_ = -1;
		resetState_ = true;
		sessionConfigured_ = false;
	}
}

void wsjtx_lib::endDecodeSession(const std::string& sessionId)
{
	std::scoped_lock lock(decodeMutex_);
	if (sessionId.empty() || sessionId == activeSessionId_) {
		activeSessionId_.clear();
		activeSessionStage_ = 0;
		activeSessionUtc_ = -1;
		stageSymbols_ = 50;
		slotUtc_ = -1;
		nagain_ = false;
		emeDelayMs_ = 0;
		resetState_ = true;
		sessionConfigured_ = false;
	}
}

bool wsjtx_lib::pullMessage(WsjtxMessage &msg) { return messageQueue_.pull(msg); }

void wsjtx_lib::decode(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int thread)
{
	std::scoped_lock lock(decodeMutex_);
	decodeLocked(mode, audiosamples, freq, thread);
}

void wsjtx_lib::decodeLocked(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int thread)
{
	std::unique_ptr<wstjx_decode> ptr = std::make_unique<wstjx_decode>();
	ptr->setStationInfo(my_call_, my_grid_, dx_call_, dx_grid_);
	ptr->setDecodeRange(decode_low_, decode_high_, decode_tol_);
	ptr->setDecodeControls(ap_decode_, decode_depth_, tx_frequency_, qso_progress_);
	ptr->setDecodeStage({stageSymbols_, slotUtc_, resetState_, nagain_, emeDelayMs_, activeSessionId_});
	wsjtx_set_message_queue(&messageQueue_);
	ptr->decode(mode, audiosamples, freq, thread);
	if (sessionConfigured_) resetState_ = false;
	wsjtx_set_message_queue(nullptr);
}

void wsjtx_lib::decode(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int thread)
{
	std::scoped_lock lock(decodeMutex_);
	decodeLocked(mode, audiosamples, freq, thread);
}

void wsjtx_lib::decodeLocked(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int thread)
{
	std::unique_ptr<wstjx_decode> ptr = std::make_unique<wstjx_decode>();
	ptr->setStationInfo(my_call_, my_grid_, dx_call_, dx_grid_);
	ptr->setDecodeRange(decode_low_, decode_high_, decode_tol_);
	ptr->setDecodeControls(ap_decode_, decode_depth_, tx_frequency_, qso_progress_);
	ptr->setDecodeStage({stageSymbols_, slotUtc_, resetState_, nagain_, emeDelayMs_, activeSessionId_});
	wsjtx_set_message_queue(&messageQueue_);
	ptr->decode(mode, audiosamples, freq, thread);
	if (sessionConfigured_) resetState_ = false;
	wsjtx_set_message_queue(nullptr);
}

void wsjtx_lib::decodeWithOptions(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int thread,
	const WsjtxDecodeConfig& config)
{
	std::scoped_lock lock(decodeMutex_);
	my_call_ = config.my_call;
	my_grid_ = config.my_grid;
	dx_call_ = config.dx_call;
	dx_grid_ = config.dx_grid;
	decode_low_ = config.low_freq;
	decode_high_ = config.high_freq;
	decode_tol_ = config.tolerance;
	ap_decode_ = config.ap_decode;
	decode_depth_ = config.decode_depth < 1 ? 1 : config.decode_depth;
	tx_frequency_ = config.tx_frequency;
	qso_progress_ = config.qso_progress < 0 ? 0 : config.qso_progress;
	setDecodeStageLocked(config.stage);
	decodeLocked(mode, audiosamples, freq, thread);
}

void wsjtx_lib::decodeWithOptions(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int thread,
	const WsjtxDecodeConfig& config)
{
	std::scoped_lock lock(decodeMutex_);
	my_call_ = config.my_call;
	my_grid_ = config.my_grid;
	dx_call_ = config.dx_call;
	dx_grid_ = config.dx_grid;
	decode_low_ = config.low_freq;
	decode_high_ = config.high_freq;
	decode_tol_ = config.tolerance;
	ap_decode_ = config.ap_decode;
	decode_depth_ = config.decode_depth < 1 ? 1 : config.decode_depth;
	tx_frequency_ = config.tx_frequency;
	qso_progress_ = config.qso_progress < 0 ? 0 : config.qso_progress;
	setDecodeStageLocked(config.stage);
	decodeLocked(mode, audiosamples, freq, thread);
}

void wsjtx_lib::decodeWithOptionsAndPull(wsjtxMode mode, WsjTxVector &audiosamples, int freq, int thread,
	const WsjtxDecodeConfig& config, std::vector<WsjtxMessage>& messages)
{
	std::scoped_lock lock(decodeMutex_);
	my_call_ = config.my_call;
	my_grid_ = config.my_grid;
	dx_call_ = config.dx_call;
	dx_grid_ = config.dx_grid;
	decode_low_ = config.low_freq;
	decode_high_ = config.high_freq;
	decode_tol_ = config.tolerance;
	ap_decode_ = config.ap_decode;
	decode_depth_ = config.decode_depth < 1 ? 1 : config.decode_depth;
	tx_frequency_ = config.tx_frequency;
	qso_progress_ = config.qso_progress < 0 ? 0 : config.qso_progress;
	setDecodeStageLocked(config.stage);
	decodeLocked(mode, audiosamples, freq, thread);
	WsjtxMessage message;
	while (messageQueue_.pull(message)) messages.push_back(std::move(message));
}

void wsjtx_lib::decodeWithOptionsAndPull(wsjtxMode mode, IntWsjTxVector &audiosamples, int freq, int thread,
	const WsjtxDecodeConfig& config, std::vector<WsjtxMessage>& messages)
{
	std::scoped_lock lock(decodeMutex_);
	my_call_ = config.my_call;
	my_grid_ = config.my_grid;
	dx_call_ = config.dx_call;
	dx_grid_ = config.dx_grid;
	decode_low_ = config.low_freq;
	decode_high_ = config.high_freq;
	decode_tol_ = config.tolerance;
	ap_decode_ = config.ap_decode;
	decode_depth_ = config.decode_depth < 1 ? 1 : config.decode_depth;
	tx_frequency_ = config.tx_frequency;
	qso_progress_ = config.qso_progress < 0 ? 0 : config.qso_progress;
	setDecodeStageLocked(config.stage);
	decodeLocked(mode, audiosamples, freq, thread);
	WsjtxMessage message;
	while (messageQueue_.pull(message)) messages.push_back(std::move(message));
}

int wspr_decode(std::vector<std::complex<float>> &iqdat, int samples,
	decoder_options options, std::vector<struct decoder_results> &decodes, int threads);

std::vector<struct decoder_results> wsjtx_lib::wspr_decode(WsjtxIQSampleVector &iqsignal, decoder_options options)
{
	std::vector<decoder_results> results;
	::wspr_decode(iqsignal, iqsignal.size(), options, results, 4);
	return results;
}

std::vector<float> wsjtx_lib::encode(wsjtxMode mode, int frequency, std::string message, std::string &messagesend, int sampleRate)
{
	switch (mode) {
	case FT8: {
		auto ptr = std::make_unique<wsjtx_encode>();
		return ptr->encode_ft8(mode, frequency, message, messagesend, sampleRate);
	}
	case FT4: {
		auto ptr = std::make_unique<wsjtx_encode>();
		return ptr->encode_ft4(mode, frequency, message, messagesend, sampleRate);
	}
	default: return {};
	}
}
