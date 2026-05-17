#include "wsjtx_encode.h"
#include "DataBuffer.h"
#include "commons.h"
#include "fortran_interface.h"
#include <algorithm>
#include <cstring>

std::vector<float> wsjtx_encode::encode_ft8(wsjtxMode mode, int frequency, std::string message, std::string &msgsent, int sampleRate)
{
	std::vector<float> signal;
	
	int i3 = 0;
	int n3 = 0;
	char ft8msgbits[77];

	std::memset(msg, 0, 38);
	std::memset(sendmsg, 0, 38);
	std::copy_n(message.c_str(), std::min<size_t>(message.size(), 37), msg);
	genft8_(msg, &i3, &n3, sendmsg, const_cast<char *>(ft8msgbits),
			const_cast<int *>(itone), 37, 37);
	sendmsg[37] = '\0';
	msgsent = std::string(sendmsg);

	int nsym = FT8_NN;
	printf("FSK tones: ");
	for (int j = 0; j < nsym; ++j)
	{
		printf("%d", itone[j]);
	}
	printf("\n");

	nsps = (sampleRate / 12000) * 1920;
	fsample = static_cast<float>(sampleRate);
	bt = 2.0;
	icmplx = 0;
	nwave = nsym * nsps;
	f0 = frequency;

	signal.clear();
	signal.resize(nwave);
	gen_ft8wave_(const_cast<int *>(itone), &nsym, &nsps, &bt, &fsample, &f0, signal.data(),
				 signal.data(), &icmplx, &nwave);
	printf("ft8 frequency %4.0f number of tones %d, samplerate %6.0f no samples %d\n", f0, nsym, fsample, nwave);
	//save_wav(signal.data(), signal.size(), FT8_SAMPLERATE, "./wave.wav");
	return signal;
}

//void genft4_(char *msg, int *ichk, char *msgsent, char ft4msgbits[], int itone[],
//			 fortran_charlen_t, fortran_charlen_t);

std::vector<float> wsjtx_encode::encode_ft4(wsjtxMode mode, int frequency, std::string message, std::string &msgsent, int sampleRate)
{
	int ichk = 0;
	char ft4msgbits[101];
	std::vector<float> signal;

	std::memset(msg, 0, 38);
	std::memset(sendmsg, 0, 38);
	std::copy_n(message.c_str(), std::min<size_t>(message.size(), 37), msg);
	genft4_(msg, &ichk, sendmsg, const_cast<char *>(ft4msgbits), const_cast<int *>(itone), 37, 37);
	sendmsg[37] = '\0';
	msgsent = std::string(sendmsg);
	
	int nsym = 103;
	int nsps = (sampleRate / 12000) * 576;
	float fsample = static_cast<float>(sampleRate);
	float f0 = frequency; //ui->TxFreqSpinBox->value() - m_XIT;
	int nwave = (nsym + 2) * nsps;
	int icmplx = 0;
	
	printf("FSK tones: ");
	for (int j = 0; j < nsym; ++j)
	{
		printf("%d", itone[j]);
	}
	printf("\n");

	signal.clear();
	signal.resize(nwave);
	gen_ft4wave_(const_cast<int *>(itone), &nsym, &nsps, &fsample, &f0, signal.data(),
				 signal.data(), &icmplx, &nwave);
	printf("ft4 frequency %d number of tones %d, samplerate %6.0f no samples %d\n", frequency, nsym, fsample, nwave);
	return signal;
}

std::vector<float> wsjtx_encode::encode_wspr(wsjtxMode mode, int frequency, std::string message, std::string &msgsent)
{
	std::vector<float> signal;

	return signal;
}
